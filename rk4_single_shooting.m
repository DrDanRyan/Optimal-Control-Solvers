function soln = rk4_single_shooting(prob, x0, tspan, uStride, varargin)

% -----------------------------------
% Initialize and setup options
% -----------------------------------
nSTATES = length(x0) + 1;
nSTEPS = length(tspan) - 1;
STATE_SHAPE = [nSTATES, nSTEPS+1];
h = tspan(2) - tspan(1); %TODO raise exception if any(diff(tspan) != 0)
tHalf = tspan(1:end-1) + h/2;

nCONTROLS = size(prob.ControlBounds, 1);
uspan = tspan(1:uStride:end); % TODO raise exception is mod(nSTEPS, uStride) !=0
nCONTROL_PTS = length(uspan);
CONTROL_SHAPE = [nCONTROLS, nCONTROL_PTS];

nodeValues = zeros(1, nCONTROL_PTS);
nodeValues(:, 2:2:end) = 1;
sawFunc = griddedInterpolant(uspan, nodeValues, 'linear');
sawLeft = sawFunc(tspan(1:end-1));
sawMid = sawFunc(tHalf);
sawRight = sawFunc(tspan(2:end));
clear nodeValues sawFunc

T0 = tspan(1);
TF = tspan(end);

if isfield(prob, 'MinMax')
   MinMax = prob.MinMax;
else
   MinMax = 'Min';
end

% u0 defaults to max(0, uMin)
u0Default = max(zeros(nCONTROLS, 1), prob.ControlBounds(:,1))*...
            ones(1, nCONTROL_PTS) + 1e-4; 

p = inputParser;
p.addParamValue('u0', u0Default);
p.addParamValue('TolX', 1e-6);
p.addParamValue('TolFun', 1e-4);
p.addParamValue('Algorithm', 'sqp');
p.addParamValue('Reporting', true);
p.addParamValue('DerivativeCheck', 'off');
parse(p, varargin{:});

u0 = p.Results.u0;
if isa(u0, 'function_handle') %if u0 is a function, eval at tspan to get vector
   u0 = u0(uspan);
end

% Turn control values into column vector for fmincon
v0 = reshape(u0, [], 1);

TolX = p.Results.TolX;
TolFun = p.Results.TolFun;
Algorithm = p.Results.Algorithm;
DerivativeCheck = p.Results.DerivativeCheck;

if p.Results.Reporting
   plotfun = @plot_func;
   iter_detail = 'iter-detailed';
else
   plotfun = [];
   iter_detail = 'off';
end


% Build nlp options structures
nlpOptions = optimoptions(@fmincon, 'Algorithm', Algorithm, ...
                                    'GradObj', 'on', ...
                                    'DerivativeCheck', DerivativeCheck, ...
                                    'TolX', TolX, ...
                                    'TolFun', TolFun, ...
                                    'Display', iter_detail, ...
                                    'PlotFcn', plotfun);


% -----------------------------------
% Main execution
% -----------------------------------

[Lb, Ub] = build_optimization_bounds();
[vOpt, soln.J] = fmincon(@nlpObjective, v0, [], [], [], [], ...
                           Lb, Ub, [], nlpOptions);
                        
if strcmp(MinMax, 'Max')
  soln.J = -soln.J;
end

[uValsOpt, uHalfValsOpt] = compute_control_values(vOpt);
xOpt = compute_states(uValsOpt, uHalfValsOpt);
lamOpt = compute_adjoints(xOpt, uValsOpt, uHalfValsOpt);

soln.u = vectorInterpolant(uspan, reshape(vOpt, CONTROL_SHAPE), 'linear');
soln.x = vectorInterpolant(tspan, xOpt(:,:,1), 'pchip');
soln.lam = vectorInterpolant(tspan, lamOpt, 'pchip');

% -----------------------------------
% Auxillary functions
% -----------------------------------
   
   function [J, dJdu] = nlpObjective(v)
      [uVals, uHalfVals] = compute_control_values(v);
      [x, J] = compute_states(uVals, uHalfVals);
      [~, dJdu] = compute_adjoints(x, uVals, uHalfVals);    
   end


   function [x, J] = compute_states(uVals, uHalfVals)
      x = nan([STATE_SHAPE, 4]); % store x, xk1, xk2, and xk3 at each time pt
      x(:,1,1) = [x0; 0];

      for i = 1:nSTEPS
         % Perform single RK-4 Step
         % F1, F2, F3, F4 refer to the RK approx values of the state rhs
         
         F1 = prob.F(tspan(i), x(:,i,1), uVals(:,i));
         x(:,i,2) = x(:,i,1) + h/2*F1;

         F2 = prob.F(tHalf(i), x(:,i,2), uHalfVals(:,i));
         x(:,i,3) = x(:,i,1) + h/2*F2;
         
         F3 = prob.F(tHalf(i), x(:,i,3), uHalfVals(:,i));
         x(:,i,4) = x(:,i,1) + h*F3;
         
         F4 = prob.F(tspan(i+1), x(:,i,4), uVals(:,i+1));
         
         x(:,i+1,1) = x(:,i,1) + h/6*(F1 + 2*F2 + 2*F3 + F4);         
      end         
      
      J = x(end,end,1);
   end


   function [lam, dJdu] = compute_adjoints(x, uVals, uHalfVals)
      lam = zeros(STATE_SHAPE);      
      lam(end,end) = 1;
      
      dJdk = nan(nSTATES, nSTEPS, 4);

      for i = nSTEPS:-1:1
         dJdk(:,i,4) = h/6*lam(:,i+1);
         dJdx3 = prob.dFdx_times_vec(tspan(i+1), x(:,i,4), uVals(:,i+1), ...
                                     dJdk(:,i,4));
                                  
         dJdk(:,i,3) = h/3*lam(:,i+1) + h*dJdx3;
         dJdx2 = prob.dFdx_times_vec(tHalf(i), x(:,i,3), uHalfVals(:,i), ...
                                     dJdk(:,i,3));
                                  
         dJdk(:,i,2) = h/3*lam(:,i+1) + h/2*dJdx2;
         dJdx1 = prob.dFdx_times_vec(tHalf(i), x(:,i,2), uHalfVals(:,i), ...
                                     dJdk(:,i,2));
                                  
         dJdk(:,i,1) = h/6*lam(:,i+1) + h/2*dJdx1;
         lam(:,i) = lam(:,i+1) + dJdx1 + dJdx2 + dJdx3 + ...
               prob.dFdx_times_vec(tspan(i), x(:,i,1), uVals(:,i), dJdk(:,i,1));
      end  
      
      if nargout > 1 % Compute dJdu
         dJdu = zeros(CONTROL_SHAPE);
         
         dJduL = prob.dFdu_times_vec(tspan(1:end-1), x(:,1:end-1,1), ...
                                       uVals(:,1:end-1), dJdk(:,:,1));
         
         dJduM = prob.dFdu_times_vec(tHalf, x(:,1:end-1,2), ...
                                       uHalfVals, dJdk(:,:,2)) ...
                 + prob.dFdu_times_vec(tHalf, x(:,1:end-1,3), ...
                                       uHalfVals, dJdk(:,:,3));
                                    
         dJduR = prob.dFdu_times_vec(tspan(2:end), x(:,1:end-1,4), ...
                                       uVals(:,2:end), dJdk(:,:,4));         
                                    
         % Compute for even index of uspan (use sawFunc)
         intervalLoss = bsxfun(@times, dJduL, sawLeft) ...
                        + bsxfun(@times, dJduM, sawMid) ...
                        + bsxfun(@times, dJduR, sawRight);
         
         % Pad ending if final control point is even index
         if mod(nCONTROL_PTS, 2) == 0
            dJdu(:,2:2:end) = squeeze(sum(reshape([intervalLoss, ...
               zeros(nCONTROLS, uStride)], nCONTROLS, 2*uStride, []), 2));
         else
            dJdu(:,2:2:end) = squeeze(sum(reshape(intervalLoss, nCONTROLS, 2*uStride, []), 2));
         end
         
         % Compute for odd index of uspan (1 - sawFunc)
         % pad loss in front since half of first saw wave is missing; transform
         % to (1 - sawFunc)*dJdu
         intervalLoss = [zeros(nCONTROLS, uStride), ...
                         dJduL + dJduM + dJduR - intervalLoss];

         % Pad ending if final control point is odd index
         if mod(nCONTROL_PTS, 2) == 1
            dJdu(:,1:2:end) = squeeze(sum(reshape([intervalLoss, ...
               zeros(nCONTROLS, uStride)], nCONTROLS, 2*uStride, []), 2));
         else
            dJdu(:,1:2:end) = squeeze(sum(reshape(intervalLoss, nCONTROLS, 2*uStride, []), 2));
         end
         
         dJdu = reshape(dJdu, [], 1);
      end
   end

   
   function [uVals, uHalfVals] = compute_control_values(v)
      u = reshape(v, CONTROL_SHAPE);      
      uFunc = vectorInterpolant(uspan, u, 'linear');
      uVals = uFunc(tspan);
      uHalfVals = uFunc(tHalf);
   end


   function [Lb, Ub] = build_optimization_bounds()
      Lb = prob.ControlBounds(:,1)*ones(1, nCONTROL_PTS);
      Ub = prob.ControlBounds(:,2)*ones(1, nCONTROL_PTS);
      Lb = reshape(Lb, [], 1);
      Ub = reshape(Ub, [], 1);
   end


   function stop = plot_func(v, optimValues, state)
      % This provides a graphical view of the progress of the NLP solver.
      % It also allows the user to stop the solver and recover the data
      % from the last iteration (as opposed to ctrl-C which will terminate
      % without returning any information)
      
      stop = 0;
      u = reshape(v, CONTROL_SHAPE);
      if strcmp(MinMax, 'Max')
         objValue = -optimValues.fval;
      else
         objValue = optimValues.fval;
      end
      
      switch state
         case 'iter'
            if optimValues.iteration == 0
               title(sprintf('Objective Value: %1.6e', objValue))
               set(gca, 'XLim', [T0 TF], 'YLim', ...
                  [min(prob.ControlBounds(:,1)), max(prob.ControlBounds(:,2))]);
               graphHandles = plot(uspan, u);
               set(graphHandles, 'Tag', 'handlesTag');
            else
               graphHandles = findobj(get(gca,'Children'),'Tag','handlesTag');
               graphHandles = graphHandles(end:-1:1); %Handle order gets reversed using findobj
               title(sprintf('Objective Value: %1.6e', objValue));
               for idx = 1:nCONTROLS
                  set(graphHandles(idx), 'YData', u(idx, :));
               end
            end
      end
   end


end

