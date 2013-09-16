function soln = single_shooting(prob, x0, controlPtArray, varargin)
% Solves the OC problem by discretizing the control as a piecewise linear 
% function and using single-shooting to compute gradient information for 
% a SQP nonlinear programming solver. 
%
% The field MinMax is optional and should have the value 'Min' or 'Max'
% depending on which type of optimization is needed.

% --------------------------------------------------------------------------------                   
% Setup constants and process optional arguments
% --------------------------------------------------------------------------------
nSTATES = length(x0);
T0 = controlPtArray(1);
TF = controlPtArray(end);
nCONTROLS = size(prob.ControlBounds, 1);
nCONTROL_PTS = length(controlPtArray);

if isfield(prob, 'MinMax')
   MinMax = prob.MinMax;
else
   MinMax = 'Min';
end


u0Default = min(prob.ControlBounds(:,2), ... % sets u0 to be zero if within ControlBounds
                max(prob.ControlBounds(:,1), zeros(nCONTROLS, 1))); % and uMin or uMax otherwise
u0Default = u0Default*ones(1, nCONTROL_PTS);

missingIC = isnan(x0); % logical vector indicating which state IC's are missing
nMISSING_IC = sum(missingIC);
ICBoundsDefault = zeros(nMISSING_IC, 2);
ICBoundsDefault(:, 1) = -Inf;
ICBoundsDefault(:, 2) = Inf;
ICGuessDefault = zeros(nMISSING_IC, 1);

p = inputParser;
p.addParamValue('u0', u0Default);
p.addParamValue('RelTol', 1e-10);
p.addParamValue('AbsTol', 1e-12);
p.addParamValue('TolX', 1e-6);
p.addParamValue('TolFun', 1e-4);
p.addParamValue('Algorithm', 'sqp');
p.addParamValue('Reporting', true);
p.addParamValue('ICBounds', ICBoundsDefault);
p.addParamValue('ICGuess', ICGuessDefault);
parse(p, varargin{:});

u0 = p.Results.u0;
if isa(u0, 'function_handle') %if u0 is a function, eval at controlPts to get vector
   u0 = u0(controlPtArray);
end
u0 = reshape(u0, [], 1);

RelTol = p.Results.RelTol;
AbsTol = p.Results.AbsTol;
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

ICBounds = p.Results.ICBounds;
v0 = [u0; p.Results.ICGuess]; % append missing state IC's to optimization vector initial guess

% Build nlp options structures
nlpOptions = optimset('Algorithm', Algorithm, ...
                      'GradObj', 'on', ...
                      'TolX', TolX, ...
                      'TolFun', TolFun, ...
                      'Display', iter_detail, ...
                      'PlotFcn', plotfun);
                   
% --------------------------------------------------------------------------------                   
% Main program
% --------------------------------------------------------------------------------
[Lb, Ub] = build_optimization_bounds();

sawFunc = build_saw_function(); % The saw function is used in the computation of 
                                % the objective function gradient.
                                
[vOpt, soln.J] = fmincon(@nlpObjective, v0, [], [], [], [], ...
                       Lb, Ub, [], nlpOptions); 

if strcmp(MinMax, 'Max')
  soln.J = -soln.J;
end

x0(missingIC) = update_IC(vOpt);
soln.u = build_control(vOpt);
[soln.x, soln.lam] = compute_x_lam(prob, x0, [T0, TF], soln.u, RelTol, AbsTol);

% --------------------------------------------------------------------------------                   
% Subfunctions
% --------------------------------------------------------------------------------

   function [objValue, dHduInt, dHduSawInt, dHdx0] = compute_obj_and_sens_integrals(v)
      % Solves ode system augmented with additional integral terms:
      % objective on the forward solve, dHdu, dHdu*sawFunc terms on the backwards solve.
      % Returns J and integral terms used for calculating dJ/dv
      
      x0(missingIC) = update_IC(v);
      u = build_control(v);
      forwardRHS_NoU = @(t, y) forwardRHS(t, y, u(t));
      [tout, yout] = odevr7(forwardRHS_NoU, [T0, TF], [x0; 0], RelTol, AbsTol);
      objValue = yout(end, end);
      x = vectorInterpolant(tout', yout(:, 1:nSTATES)', 'pchip');
      
      backwardInitCond = zeros(nSTATES+2*nCONTROLS, 1);
      backwardRHS_NoXU = @(t, y) backwardRHS(t, x(t), y, u(t));
      [tout, yout] = odevr7(backwardRHS_NoXU, [TF, T0], backwardInitCond, ...
                               RelTol, AbsTol);
                            
      dHduIntFunc = vectorInterpolant(flip(tout', 2), ...
                           flip(yout(:, nSTATES+1:nSTATES+nCONTROLS)', 2), 'pchip');
      dHduSawIntFunc = vectorInterpolant(flip(tout', 2), ...
                           flip(yout(:, nSTATES+nCONTROLS+1:end)', 2), 'pchip');
      
      dHduInt = dHduIntFunc(controlPtArray);
      dHduSawInt = dHduSawIntFunc(controlPtArray);
      dHdx0 = yout(end, [missingIC; false; false]); % dHdx0 equal adjoint values at t=0 for missing state IC's
   end

   function [objValue, gradValue] = nlpObjective(v)
      % Computes objective functional value and gradient of objective
      % functional.
      
      [objValue, dHduInt, dHduSawInt, dHdx0] = compute_obj_and_sens_integrals(v);

      gradValue = zeros(nCONTROLS*nCONTROL_PTS + nMISSING_IC, 1);
      for idx = 1:nCONTROL_PTS
         left_idx = max(1, idx-1);
         right_idx = min(nCONTROL_PTS, idx+1);
         gradStartIdx = nCONTROLS*(idx-1)+1;
         gradStopIdx = nCONTROLS*idx;
         if mod(idx,2) % idx is odd
            gradValue(gradStartIdx:gradStopIdx) = dHduInt(:, right_idx) - dHduInt(:, left_idx) ...
                                                  - (dHduSawInt(:, right_idx) - dHduSawInt(:, left_idx));
         else % idx is even
            gradValue(gradStartIdx:gradStopIdx) = dHduSawInt(:, right_idx) - dHduSawInt(:, left_idx);
         end
      end
      
      if nMISSING_IC > 0
         gradValue(nCONTROLS*nCONTROL_PTS+1:end) = dHdx0;
      end
      
      if strcmp(MinMax, 'Max')
         objValue = -objValue;
         gradValue = -gradValue;
      end
   end
   
   function value = forwardRHS(t, y, u)
      % The state equation augmented with the objective integral
      tLENGTH = length(t);
      value = zeros(nSTATES+1, tLENGTH);
      x = y(1:nSTATES, :);
      value(1:nSTATES, :) = prob.stateRHS(t, x, u);      
      value(end, :) = prob.objective(t, x, u);
   end

   function value = backwardRHS(t, x, y, u)
      % The adjoint equation augmented with an integral of dHdu and
      % dHdu*sawFunc, both of which are used to compute gradient of
      % objective function.
      tLENGTH = size(t, 2);
      value = zeros(nSTATES+2*nCONTROLS, tLENGTH);
      lam = y(1:nSTATES, :);
      
      value(1:nSTATES, :) = prob.adjointRHS(t, x, lam, u);
      dHduValue = prob.dHdu(t, x, lam, u);
      value(nSTATES+1:end, :) = [dHduValue; dHduValue.*sawFunc(t)];
   end

   function [Lb, Ub] = build_optimization_bounds()
      Lb = prob.ControlBounds(:,1)*ones(1, nCONTROL_PTS);
      Ub = prob.ControlBounds(:,2)*ones(1, nCONTROL_PTS);
      Lb = reshape(Lb, [], 1);
      Ub = reshape(Ub, [], 1);
      Lb = [Lb; ICBounds(:,1)];
      Ub = [Ub; ICBounds(:,2)];
   end

   function sawFunc = build_saw_function()
      nodeValues = zeros(nCONTROLS, nCONTROL_PTS);
      nodeValues(:, 2:2:end) = 1;
      sawFunc = vectorInterpolant(controlPtArray, nodeValues, 'linear');
   end

   function control_func = build_control(v)
      v = reshape(v(1:nCONTROLS*nCONTROL_PTS), nCONTROLS, []);
      control_func = vectorInterpolant(controlPtArray, v, 'linear');
   end

   function missingVals = update_IC(v)
      if nMISSING_IC == 0
         missingVals = [];
      else
         missingVals = v(nCONTROLS*nCONTROL_PTS+1:end);
      end
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

end