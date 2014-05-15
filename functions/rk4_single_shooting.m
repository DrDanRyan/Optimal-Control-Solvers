function soln = rk4_single_shooting(prob, x0, tspan, nCONTROL_PTS, varargin)

% -----------------------------------
% Initialize and setup options
% -----------------------------------
nSTATES = length(x0) + 1;
nSTEPS = length(tspan) - 1;
STATE_SHAPE = [nSTATES, nSTEPS+1];
h = diff(tspan);

t = zeros(1, 2*nSTEPS+1); % a time vector that also includes the interval midpoints
t(1:2:end) = tspan;
t(2:2:end-1) = (tspan(1:end-1) + tspan(2:end))/2;

T0 = tspan(1);
TF = tspan(end);

nCONTROLS = size(prob.ControlBounds, 1);

if isfield(prob, 'MinMax')
   MinMax = prob.MinMax;
else
   MinMax = 'Min';
end

p = inputParser;
p.addParamValue('TolX', 1e-6);
p.addParamValue('TolFun', 1e-4);
p.addParamValue('Algorithm', 'sqp');
p.addParamValue('Reporting', true);
p.addParamValue('DerivativeCheck', 'off');
p.addParamValue('ControlType', 'linear', ...
                  @(x) any(strcmp(x, {'linear', 'Chebyshev'})));
parse(p, varargin{:});

% Turn control values into column vector for fmincon
TolX = p.Results.TolX;
TolFun = p.Results.TolFun;
Algorithm = p.Results.Algorithm;
DerivativeCheck = p.Results.DerivativeCheck;
ControlType = p.Results.ControlType;

switch ControlType
   case 'linear'
      control = PWLinearControl(t, nCONTROL_PTS, nCONTROLS);
   case 'Chebyshev'
      control = ChebyshevControl(t, nCONTROL_PTS, nCONTROLS);
end

if p.Results.Reporting
   plotfun = @plot_func;
   iter_detail = 'iter-detailed';
else
   plotfun = [];
   iter_detail = 'off';
end

% If approximating infinite time solution, canonical equilibrium and stable
% eigenvectors. Then precompute Vlam*Vx^{-1} (see notes).
if isInfiniteApprox
   
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

v0 = control.compute_initial_v(prob.ControlBounds);

if ismethod(control, 'compute_nlp_bounds')
   [Lb, Ub] = control.compute_nlp_bounds(prob.ControlBounds);
else
   Lb = [];
   Ub = [];
end

if ismethod(control, 'compute_nonlcon')
   nonlcon = @control.compute_nonlcon;
else
   nonlcon = [];
end

[vOpt, soln.J] = fmincon(@nlpObjective, v0, [], [], [], [], ...
                           Lb, Ub, nonlcon, nlpOptions);
                        
if strcmp(MinMax, 'Max')
  soln.J = -soln.J;
end

uOpt = control.compute_u(t, vOpt);
xOpt = compute_states(uOpt);
lamOpt = compute_adjoints(xOpt, uOpt);

soln.u = control.compute_uFunc(vOpt);
soln.x = vectorInterpolant(tspan, xOpt(:,:,1), 'pchip');
soln.lam = vectorInterpolant(tspan, lamOpt, 'pchip');

% -----------------------------------
% Auxillary functions
% -----------------------------------
   
   function [J, dJdv] = nlpObjective(v)
      u = control.compute_u(t, v);
      [x, J] = compute_states(u);
      [~, dJdv] = compute_adjoints(x, u);    
   end


   function [x, J] = compute_states(u)
      x = nan([STATE_SHAPE, 4]); % store x, xk1, xk2, and xk3 at each time pt
      x(:,1,1) = [x0; 0];

      for i = 1:nSTEPS
         % Perform single RK-4 Step
         % F1, F2, F3, F4 refer to the RK approx values of the state rhs
         
         F1 = prob.F(t(2*i-1), x(:,i,1), u(:,2*i-1));
         x(:,i,2) = x(:,i,1) + h(i)/2*F1;

         F2 = prob.F(t(2*i), x(:,i,2), u(:,2*i));
         x(:,i,3) = x(:,i,1) + h(i)/2*F2;
         
         F3 = prob.F(t(2*i), x(:,i,3), u(:,2*i));
         x(:,i,4) = x(:,i,1) + h(i)*F3;
         
         F4 = prob.F(t(2*i+1), x(:,i,4), u(:,2*i+1));
         
         x(:,i+1,1) = x(:,i,1) + h(i)/6*(F1 + 2*F2 + 2*F3 + F4);         
      end         
      
      J = x(end,end,1);
   end


   function [lam, dJdv] = compute_adjoints(x, u)
      
      lam = zeros(STATE_SHAPE);      
      lam(end,end) = 1;
            
      dJdk = nan(nSTATES, nSTEPS, 4);

      for i = nSTEPS:-1:1
         dJdk(:,i,4) = h(i)/6*lam(:,i+1);
         dJdx3 = prob.dFdx_times_vec(t(2*i+1), x(:,i,4), u(:,2*i+1), ...
                                     dJdk(:,i,4));
                                  
         dJdk(:,i,3) = h(i)/3*lam(:,i+1) + h(i)*dJdx3;
         dJdx2 = prob.dFdx_times_vec(t(2*i), x(:,i,3), u(:,2*i), ...
                                     dJdk(:,i,3));
                                  
         dJdk(:,i,2) = h(i)/3*lam(:,i+1) + h(i)/2*dJdx2;
         dJdx1 = prob.dFdx_times_vec(t(2*i), x(:,i,2), u(:,2*i), ...
                                     dJdk(:,i,2));
                                  
         dJdk(:,i,1) = h(i)/6*lam(:,i+1) + h(i)/2*dJdx1;
         lam(:,i) = lam(:,i+1) + dJdx1 + dJdx2 + dJdx3 + ...
               prob.dFdx_times_vec(t(2*i-1), x(:,i,1), u(:,2*i-1), dJdk(:,i,1));
      end  
      
      if nargout > 1
         dJdu = compute_dJdu(x, u, dJdk);
         dJdv = control.compute_dJdv(dJdu);
      end
   end

   
   function dJdu = compute_dJdu(x, u, dJdk)
         dJdu = zeros(size(u)); % ~ nCONTROLS x 2*nSTEPS+1
         
         % Left end point
         dJdu(:,1) = prob.dFdu_times_vec(t(1), x(:,1,1), u(:,1), dJdk(:,1,1));
         
         % RK step interval mid points
         dJdu(:,2:2:end-1) = ...
            prob.dFdu_times_vec(t(2:2:end-1), x(:,1:end-1,2), ...
                                u(:,2:2:end-1), dJdk(:,:,2)) ...
            + prob.dFdu_times_vec(t(2:2:end-1), x(:,1:end-1,3), ...
                                  u(:,2:2:end-1), dJdk(:,:,3));
         
         % RK step interval endpts (except for left and right of full interval
         dJdu(:,3:2:end-2) = ...
            prob.dFdu_times_vec(t(3:2:end-2), x(:,2:end-1,1), ...
                                u(:,3:2:end-2), dJdk(:,2:end,1)) ...
            + prob.dFdu_times_vec(t(3:2:end-2), x(:,1:end-2,4), ...
                                  u(:,3:2:end-2), dJdk(:,1:end-1,4));
         
         % Right end point
         dJdu(:,end) = prob.dFdu_times_vec(t(end), x(:,end-1,4), ...
                                           u(:,end), dJdk(:,end,4)); 
   end


   function stop = plot_func(v, optimValues, state)
      % This provides a graphical view of the progress of the NLP solver.
      % It also allows the user to stop the solver and recover the data
      % from the last iteration (as opposed to ctrl-C which will terminate
      % without returning any information)
      
      stop = 0;
      u = control.compute_u(tspan, v);
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

