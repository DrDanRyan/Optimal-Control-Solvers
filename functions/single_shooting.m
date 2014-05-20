function soln = single_shooting(prob, x0, tspan, nCONTROL_PTS, varargin)

% -----------------------------------
% Initialize and setup options
% -----------------------------------
T0 = tspan(1);
TF = tspan(end);

nCONTROLS = size(prob.ControlBounds, 1);

if isfield(prob, 'MinMax')
   MinMax = prob.MinMax;
else
   MinMax = 'Min';
end


% Parse varargin
p = inputParser;
p.addParamValue('TolX', 1e-5);
p.addParamValue('TolFun', 3e-4);
p.addParamValue('Algorithm', 'sqp');
p.addParamValue('Reporting', true);
p.addParamValue('DerivativeCheck', 'off');
p.addParamValue('Control', []);
p.addParamValue('Integrator', []);
p.addParamValue('u0', 0);
parse(p, varargin{:});

TolX = p.Results.TolX;
TolFun = p.Results.TolFun;
Algorithm = p.Results.Algorithm;
DerivativeCheck = p.Results.DerivativeCheck;


% Initialize integrator and control objects
if isempty(p.Results.Integrator)
   integrator = RK4Integrator(tspan);
else
   integrator = p.Results.Integrator;
end


if isempty(p.Results.Control)
   control = PWLinearControl(integrator.t, nCONTROL_PTS, nCONTROLS);
else
   control = p.Results.Control;
end


% Initialize initial control values
u0 = min(prob.ControlBounds(:,2), max(prob.ControlBounds(:,1), p.Results.u0));


% Setup plot function
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

v0 = control.compute_initial_v(u0);

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

uOpt = control.compute_u(vOpt);
xOpt = integrator.compute_states(prob, x0, uOpt);
lamOpt = integrator.compute_adjoints(prob, uOpt);

soln.u = control.compute_uFunc(vOpt);
soln.x = vectorInterpolant(tspan, xOpt(1:end-1,:), 'pchip');
soln.lam = vectorInterpolant(tspan, lamOpt(1:end-1,:), 'pchip');


% -----------------------------------
% Auxillary functions
% -----------------------------------
   
   function [J, dJdv] = nlpObjective(v)
      u = control.compute_u(v);
      [~, J] = integrator.compute_states(prob, x0, u);
      [~, dJdu] = integrator.compute_adjoints(prob, u);    
      dJdv = control.compute_dJdv(dJdu);
   end


   function stop = plot_func(v, optimValues, state)
      % This provides a graphical view of the progress of the NLP solver.
      % It also allows the user to stop the solver and recover the data
      % from the last iteration (as opposed to ctrl-C which will terminate
      % without returning any information)
      
      stop = 0;
      u = control.compute_u(v);
      if strcmp(MinMax, 'Max')
         objValue = -optimValues.fval;
      else
         objValue = optimValues.fval;
      end
      
      if strcmp(state, 'iter')
         if optimValues.iteration == 0
            title(sprintf('Objective Value: %1.6e', objValue))
            set(gca, 'XLim', [T0 TF], 'YLim', ...
               [min(prob.ControlBounds(:,1)), max(prob.ControlBounds(:,2))]);
            graphHandles = plot(integrator.t, u);
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

