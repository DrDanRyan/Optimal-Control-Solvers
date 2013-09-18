function soln = bvp_solver(prob, x0, tspan, varargin)
% Solves the optimal control problem using one of Matlab's native bvp
% solvers (default is bvp5c). Can accept initial guesses for either 
% y = [x; lam] OR the control function. Either way, the guess should be in
% the form of a function handle.

% --------------------------------------------------------------------------------                   
% Setup constants and process options structure
% --------------------------------------------------------------------------------
T0 = tspan(1);
TF = tspan(end);
nSTATES = length(x0); % this is revised below in the case of user-supplied boundary conditions

% Assign default values to optional arguments
bvpRelTol = 1e-6;
bvpAbsTol = 1e-6;
RelTol = 1e-8;
AbsTol = 1e-8;
MAX_MESH_SIZE = 1000;
nINTERP_PTS = 1001; 
Solver = @bvp5c;
InitMesh = linspace(T0, TF, 11);


% Overwrite defaults if user supplied value
if size(varargin) > 0
   options = varargin{1};
else
   options = struct();
end
   
if isfield(options, 'bvpRelTol')
   bvpRelTol = options.bvpRelTol;
end

if isfield(options, 'bvpAbsTol')
   bvpAbsTol = options.bvpAbsTol;
end

if isfield(options, 'RelTol')
   RelTol = options.RelTol;
end

if isfield(options, 'AbsTol')
   AbsTol = options.AbsTol;
end

if isfield(options, 'MAX_MESH_SIZE')
   MAX_MESH_SIZE = options.MAX_MESH_SIZE;
end

if isfield(options, 'nINTERP_PTS')
   nINTERP_PTS = options.nINTERP_PTS;
end

if isfield(options, 'bvpSolver')
   Solver = options.bvpSolver;
end

if isfield(options, 'InitMesh')
   InitMesh = options.InitMesh;
end


% Define boundary condition functions
bcFunc = @(yL, yR) [yL(1:end/2) - x0; yR(end/2+1:end)];

yLJac = zeros(2*nSTATES);
yRJac = zeros(2*nSTATES);
yLJac(1:nSTATES, 1:nSTATES) = eye(nSTATES);
yRJac(nSTATES+1:end, nSTATES+1:end) = eye(nSTATES);

% Overwrite boundary conditions if provided in problem structure
if isfield(prob, 'bcFunc')
   bcFunc = prob.bcFunc;
end

if isfield(prob, 'yLJac')
   yLJac = prob.yLJac;
end

if isfield(prob, 'yRJac')
   yRJac = prob.yRJac;
end

nSTATES = size(yLJac, 1)/2; %used to set correct number of states when custom bc's are used
y0vec = zeros(2*nSTATES, 1);
y0vec(1:length(x0)) = x0; %if there are states without initial conditions, they are set to zero
y0 = bvpinit(InitMesh, y0vec);

if isfield(options, 'y0')
   yGuess = options.y0;   
   y0 = bvpinit(InitMesh, yGuess);
elseif isfield(options, 'u0')
   [xGuess, lamGuess] = compute_x_lam(prob, x0, tspan, options.u0, RelTol, AbsTol);
   yGuess = @(t) [xGuess(t); lamGuess(t)];
   y0 = bvpinit(InitMesh, yGuess);
end
interpPts = linspace(T0, TF, nINTERP_PTS);

% --------------------------------------------------------------------------------                   
% Main Program
% --------------------------------------------------------------------------------



% Define optimality system
   function value = optRHS(t, y)
      x = y(1:end/2, :); 
      lam = y(end/2+1:end, :); 
      u = prob.ControlChar(t, x, lam);
      value = [prob.stateRHS(t, x, u); prob.adjointRHS(t, x, lam, u)];
   end

% Build bvp options structure
bvpOptions = bvpset('RelTol', bvpRelTol, 'AbsTol', bvpAbsTol, 'NMax', MAX_MESH_SIZE, ...
                       'BCJacobian', {yLJac, yRJac}, 'Vectorized', 'on');

if isfield(prob, 'optJac')
   bvpOptions = bvpset(bvpOptions, 'FJacobian', prob.optJac);
end

% Solve the BVP problem
y = Solver(@optRHS, bcFunc, y0, bvpOptions);

% Make a optimal control function using interpolation (so that the control
% does not directly depend on any global parameters).
xVals = deval(y, interpPts, 1:nSTATES);
lamVals = deval(y, interpPts, nSTATES+1:2*nSTATES);
soln.u = vectorInterpolant(interpPts, prob.ControlChar(interpPts, xVals, lamVals), 'pchip');

% Recompute x and lambda using odevr7 so that they are in separate
% structures (soln should have same interface as other solvers). Also
% computes objective value in the process.
x0 = xVals(:, 1); % needed in case of user supplied bc's
[soln.x, soln.lam, soln.J] = compute_x_lam_J(prob, x0, tspan, soln.u, RelTol, AbsTol);

end