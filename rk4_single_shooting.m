function soln = rk4_single_shooting(prob, x0, tspan, varargin)

% -----------------------------------
% Initialize and setup options
% -----------------------------------
nSTATES = length(x0);
nSTEPS = length(tspan) - 1;
h = tspan(2) - tspan(1); %TODO raise exception if ~all(diff(tspan) == 0)
nCONTROLS = size(prob.ControlBounds, 1);

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

soln.u = build_control(vOpt);
[soln.x, soln.lam] = compute_x_lam(prob, x0, [T0, TF], soln.u);


% -----------------------------------
% Auxillary functions
% -----------------------------------
   
   function [J, dJdu] = nlpObjective(v)
      
   end


   function [Lb, Ub] = build_optimization_bounds()
      Lb = prob.ControlBounds(:,1)*ones(1, nCONTROL_PTS);
      Ub = prob.ControlBounds(:,2)*ones(1, nCONTROL_PTS);
      Lb = reshape(Lb, [], 1);
      Ub = reshape(Ub, [], 1);
   end


   function control_func = build_control(v)
      v = reshape(v(1:nCONTROLS*nCONTROL_PTS), nCONTROLS, []);
      control_func = vectorInterpolant(controlPtArray, v, 'linear');
   end


   function stop = plot_func(v, optimValues, state)
      % This provides a graphical view of the progress of the NLP solver.
      % It also allows the user to stop the solver and recover the data
      % from the last iteration (as opposed to ctrl-C which will terminate
      % without returning any information)
      
      stop = 0;
      uValues = reshape(v(1:nCONTROLS*nCONTROL_PTS), nCONTROLS, []);
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
               graphHandles = plot(controlPtArray, uValues);
               set(graphHandles, 'Tag', 'handlesTag');
            else
               graphHandles = findobj(get(gca,'Children'),'Tag','handlesTag');
               graphHandles = graphHandles(end:-1:1); %Handle order gets reversed using findobj
               title(sprintf('Objective Value: %1.6e', objValue));
               for idx = 1:nCONTROLS
                  set(graphHandles(idx), 'YData', uValues(idx, :));
               end
            end
      end
   end


   function midPts = compute_midpoints(vec)
      midPts = .5*(vec(:,1:end-1) + vec(:,2:end));
   end


   function [x, J, xK] = compute_states(u)
      x = zeros(nStates, nSteps+1);
      x(:,1) = x0;
      xK = zeros(nStates, nSteps, 3);
      cumObj = zeros(1, nSteps+1);

      h = h;
      tHalf = compute_midpoints(t);
      uHalf = compute_midpoints(u);

      for i = 1:nSteps
         y1 = f(t(i), x(:,i), u(:,i));
         xK(:,i,1) = x(:,i) + .5*h*y1;

         y2 = f(tHalf(i), xK(:,i,1), uHalf(:,i));
         xK(:,i,2) = x(:,i) + .5*h*y2;

         y3 = f(tHalf(i), xK(:,i,2), uHalf(:,i));
         xK(:,i,3) = x(:,i) + h*y3;

         y4 = f(t(i+1), xK(:,i,3), u(:,i+1));

         x(:,i+1) = x(:,i) + h*(y1 + 2*y2 + 2*y3 + y4)/6;

         y1 = g(t(i), x(:,i), u(:,i));
         y2 = g(tHalf(i), xK(:,i,1), uHalf(:,i));
         y3 = g(tHalf(i), xK(:,i,2), uHalf(:,i));
         y4 = g(t(i+1), xK(:,i,3), u(:,i+1));

         cumObj(i+1) = cumObj(i) + h*(y1 + 2*y2 + 2*y3 + y4)/6;
      end
      J = cumObj(end);         
   end


   function lam = compute_adjoints(obj, x, xK, u)
      lam = zeros(nStates, nSteps+1);
      t = tspan;
      h = h;
      tHalf = compute_midpoints(t);
      uHalf = compute_midpoints(u);

      for i = nSteps:-1:1
         y1 = dfdx(t(i), x(:,i), u(:,i));
         y2 = dfdx(tHalf(i), xK(:,i,1), uHalf(:,i))*...
            (eye(nStates) + .5*h*y1);  % TODO check if matrix order is correct in multiplication
         y3 = dfdx(tHalf(i), xK(:,i,2), uHalf(:,i))*...
            (eye(nStates) + .5*h*y2);
         y4 = dfdx(t(i+1), xK(:,i,3), u(:,i+1))*...
            (eye(nStates) + h*y3);


      end
   end

end

