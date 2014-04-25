function soln = rk4_single_shooting(prob, x0, tspan, varargin)

% -----------------------------------
% Initialize and setup options
% -----------------------------------
nSTATES = length(x0);
nSTEPS = length(tspan) - 1;
STATE_SHAPE = [nSTATES, nSTEPS+1];
h = tspan(2) - tspan(1); %TODO raise exception if ~all(diff(tspan) == 0)
nCONTROLS = size(prob.ControlBounds, 1);
CONTROL_SHAPE = [nCONTROLS, nSTEPS+1];

T0 = tspan(1);
TF = tspan(end);

if isfield(prob, 'MinMax')
   MinMax = prob.MinMax;
else
   MinMax = 'Min';
end

% u0 defaults to max(0, uMin)
u0Default = max(zeros(nCONTROLS, 1), prob.ControlBounds(:,1))*...
            ones(1, nSTEPS+1); 

p = inputParser;
p.addParamValue('u0', u0Default);
p.addParamValue('TolX', 1e-6);
p.addParamValue('TolFun', 1e-4);
p.addParamValue('Algorithm', 'sqp');
p.addParamValue('Reporting', true);
parse(p, varargin{:});

u0 = p.Results.u0;
if isa(u0, 'function_handle') %if u0 is a function, eval at tspan to get vector
   u0 = u0(tspan);
end

% Turn control values into column vector for fmincon
v0 = reshape(u0, [], 1);

TolX = p.Results.TolX;
TolFun = p.Results.TolFun;
Algorithm = p.Results.Algorithm;

if p.Results.Reporting
   plotfun = @plot_func;
   iter_detail = 'iter-detailed';
else
   plotfun = [];
   iter_detail = 'off';
end


% Build nlp options structures
nlpOptions = optimset('Algorithm', Algorithm, ...
                      'GradObj', 'on', ...
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

uOpt = reshape(vOpt, CONTROL_SHAPE);
xOpt = compute_states(uOpt);
lamOpt = compute_adjoints(xOpt, uOpt);

soln.u = vectorInterpolant(tspan, uOpt, 'linear');
soln.x = vectorInterpolant(tspan, xOpt(:,:,1), 'pchip');
soln.lam = vectorInterpolant(tspan, lamOpt, 'pchip');

% -----------------------------------
% Auxillary functions
% -----------------------------------
   
   function [J, dJdu] = nlpObjective(v)
      u = reshape(v, CONTROL_SHAPE);
      [x, J] = compute_states(u);
      [~, dJdu] = compute_adjoints(x, u);    
   end


   function [x, J] = compute_states(u)
      x = nan([STATE_SHAPE, 4]); % store x, xk1, xk2, and xk3 at each time pt
      x(:,1,1) = x0;

      tHalf = compute_midpoints(tspan);
      uHalf = compute_midpoints(u);

      for i = 1:nSTEPS
         % Perform single RK-4 Step
         % F1, F2, F3, F4 refer to the RK approx values of the state rhs
         
         F1 = prob.F(tspan(i), x(:,i,1), u(:,i));
         x(:,i,2) = x(:,i,1) + h/2*F1;

         F2 = prob.F(tHalf(i), x(:,i,2), uHalf(:,i));
         x(:,i,3) = x(:,i,1) + h/2*F2;
         
         F3 = prob.F(tHalf(i), x(:,i,3), uHalf(:,i));
         x(:,i,4) = x(:,i,1) + h*F3;
         
         F4 = prob.F(tspan(i+1), x(:,i,4), u(:,i+1));
         
         x(:,i+1,1) = x(:,i,1) + h/6*(F1 + 2*F2 + 2*F3 + F4);         
      end         
      
      J = x(end,end,1);
   end


   function [lam, dJdu] = compute_adjoints(x, u)
      
      lam = zeros(STATE_SHAPE);      
      lam(end,end) = 1;

      tHalf = compute_midpoints(tspan);
      uHalf = compute_midpoints(u);
      
      dJdk = nan(nSTATES, nSTEPS, 4);

      for i = nSTEPS:-1:1
         dJdk(:,i,4) = h/6*lam(:,i+1);
         dJdx3 = prob.dFdx_times_vec(tspan(i+1), x(:,i,4), u(:,i+1), ...
                                     dJdk(:,i,4));
                                  
         dJdk(:,i,3) = h/3*lam(:,i+1) + h*dJdx3;
         dJdx2 = prob.dFdx_times_vec(tHalf(i), x(:,i,3), uHalf(i), ...
                                     dJdk(:,i,3));
                                  
         dJdk(:,i,2) = h/3*lam(:,i+1) + h/2*dJdx2;
         dJdx1 = prob.dFdx_times_vec(tHalf(i), x(:,i,2), uHalf(i), ...
                                     dJdk(:,i,2));
                                  
         dJdk(:,i,1) = h/6*lam(:,i+1) + h/2*dJdx1;
         lam(:,i) = lam(:,i+1) + dJdx1 + dJdx2 + dJdx3 + ...
                   prob.dFdx_times_vec(tspan(i), x(:,i,1), u(:,i), dJdk(:,i,1));
      end  
      
      if nargout > 1 % Compute dJdu
         dJdu = zeros(CONTROL_SHAPE);
         
         leftend = prob.dFdu_times_vec(tspan(1:end-1), x(:,1:end-1,1), ...
                                       u(:,1:end-1), dJdk(:,:,1));
         
         leftmid = .5*prob.dFdu_times_vec(tHalf, x(:,1:end-1,2), ...
                                       uHalf, dJdk(:,:,2));
            
         rightmid = .5*prob.dFdu_times_vec(tHalf, x(:,1:end-1,3), ...
                                       uHalf, dJdk(:,:,3));
                                    
         rightend = prob.dFdu_times_vec(tspan(2:end), x(:,1:end-1,4), ...
                                       u(:,2:end), dJdk(:,:,4));
         
         % Compute contributions from the right of u_j
         dJdu(:,1:end-1) = leftend + leftmid + rightmid;
         
         % Compute contributions from the left of u_j
         dJdu(:,2:end) = dJdu(:,2:end) + leftmid + rightmid + rightend;
      end
   end


   function [Lb, Ub] = build_optimization_bounds()
      Lb = prob.ControlBounds(:,1)*ones(1, nSTEPS+1);
      Ub = prob.ControlBounds(:,2)*ones(1, nSTEPS+1);
      Lb = reshape(Lb, [], 1);
      Ub = reshape(Ub, [], 1);
   end


   function midPts = compute_midpoints(vec)
      midPts = .5*(vec(:,1:end-1) + vec(:,2:end));
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
               graphHandles = plot(tspan, u);
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

