function soln = fb_sweep(prob, x0, tspan, varargin)
% Attempts to solve an optimal control problem by forward backwards sweep.
% May fail to converge in which case soln will be an empty struct.
%
% If provided, u0 should be either a vector of equally spaced control
% values or a function of a single variable (time).


% --------------------------------------------------------------------------------                   
% Setup constants and process options structure
% --------------------------------------------------------------------------------
T0 = tspan(1);
TF = tspan(end);

% Set default values for optional arguments
uRelTol = 1e-7; % Relative tolerance for change in control for a single sweep
uAbsTol = 1e-7; % Absolute tolerance for change in control for a single sweep 
RelTol = 5e-14;
AbsTol = 5e-14;
nSWEEPS = 50;
nERROR_PTS = 1001; % Number of points where change in control is measured
nINTERP_PTS = 1001;
u0 = @(t) prob.ControlBounds(:,1)*ones(1, length(t)); % set u0 to lower bound(s) of control(s)

% Overwrite default optional values if found in options structure
if size(varargin) > 0
   options = varargin{1};
else
   options = struct();
end

if isfield(options, 'uRelTol')
   uRelTol = options.uRelTol;
end

if isfield(options, 'uAbsTol')
   uAbsTol = options.uAbsTol;
end

if isfield(options, 'RelTol')
   RelTol = options.RelTol;
end

if isfield(options, 'AbsTol')
   AbsTol = options.AbsTol;
end

if isfield(options, 'nSWEEPS')
   nSWEEPS = options.nSWEEPS;
end

if isfield(options, 'nERROR_PTS')
   nERROR_PTS = options.nERROR_PTS;
end

if isfield(options, 'nINTERP_PTS')
   nINTERP_PTS = options.nINTERP_PTS;
end

if isfield(options, 'u0')
   u0 = options.u0;
   if isa(u0, 'numeric') % if u0 a vector/array then turn u0 into a function
   % assume u0 values are evenly spaced between T0 and TF and build interpolant function
      time = linspace(T0, TF, size(u0, 2));
      u0 = vectorInterpolant(time, u0, 'pchip');
   end
end

errorPts = linspace(T0, TF, nERROR_PTS);
interpPts = linspace(T0, TF, nINTERP_PTS);
   

% --------------------------------------------------------------------------------                   
% Main Program
% --------------------------------------------------------------------------------
u = u0;
soln = struct(); % will remain empty if convergence fails

for sweepIdx = 1:nSWEEPS
   uNew = sweep(u);
   if check_convergence(uNew, u)
      [soln.x, soln.lam, soln.u, soln.J] = final_sweep(u);
      break;
   else
      u = uNew;
   end
end


% --------------------------------------------------------------------------------                   
% Subfunctions
% --------------------------------------------------------------------------------

   function uNew = sweep(u)
      [x, lam] = compute_x_lam(prob, x0, [T0, TF], u, RelTol, AbsTol);
      uNew = @(t) prob.ControlChar(t, x(t), lam(t));
   end

   function flag = check_convergence(uNew, u)
      % Computes a weighted difference vecotr between the new and old
      % control functions at each point in errorPts. If the value is < 1,
      % then the difference satisfies the tolerance for the change in u.
      % The algorithm will terminate if the tolerance is satisfied at every
      % errorPt. The maximum of the computed values is printed to the
      % screen.
      
      weightedVec = abs(uNew(errorPts) - u(errorPts))./(uRelTol*abs(u(errorPts)) + uAbsTol);
      maxChange = max(max(weightedVec));
      fprintf('Iteration: %d \t Normalized change in u: %f \n', sweepIdx, maxChange)
      if maxChange <= 1
         flag = true;
      else 
         flag = false;
      end
   end

   function [xOpt, lamOpt, uOpt, objValue] = final_sweep(u)
      % One last sweep that also computes the objective value.
      % uOpt is interpolated so that it does not depend on any global
      % parameters.
      
      [xOpt, lamOpt, objValue] = compute_x_lam_J(prob, x0, [T0, TF], u, RelTol, AbsTol);
      uOpt = vectorInterpolant(interpPts, prob.ControlChar(interpPts, xOpt(interpPts), ...
                                                            lamOpt(interpPts)), 'pchip');        
   end

end