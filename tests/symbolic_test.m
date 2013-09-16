function prob = symbolic_test()

nStates = 2;
nControls = 2;
x = sym('x', [1 nStates]);
lam = sym('lam', [1 nStates]);
u = sym('u', [1, nControls]);

% Objective integrand
obj(x, u) = x(1)^2 + x(2)^2 + u(1)^2 + u(2)^2;

% rhs of ODE system
f(x, u) = [x(1)*x(2) - u(1); 
           u(2)*x(2) + 3];

% Hamiltonian
H(x, u, lam) = obj + lam*f;
 
% % Adjoint ODE rhs
g = -gradient(H, x);

% dH/du
dHdu = gradient(H, u);

% Optimality condition
uCellTemp = unpack(u);
uOptCell = cell(1, nControls);
[uOptCell{:}] = solve(dHdu, uCellTemp{:});
uOpt(x, lam) = transpose(sym(uOptCell));
clear uCell uOptCell

% Convert symbolic function to matlab function handles
obj = matlabFunction(obj);
f = matlabFunction(f);
H = matlabFunction(H);
g = matlabFunction(g);
dHduTemp = matlabFunction(dHdu);
uOpt = matlabFunction(uOpt);

prob.objective = @objective;
prob.stateRHS = @stateRHS;

   function value = objective(t, x, u)
      xCell = cell(1, nStates);
      uCell = cell(1, nControls);
      for i = 1:nStates
         xCell{i} = x(i,:);
      end

      for i = 1:nControls
         uCell{i} = u(i,:);
      end

      value = obj(xCell{:}, uCell{:});
   end

   function value = stateRHS(t, x, u)
      xCell = cell(1, nStates);
      uCell = cell(1, nControls);
      for i = 1:nStates
         xCell{i} = x(i, :);
      end

      for i = 1:nControls
         uCell{i} = u(i,:);
      end

      value = f(xCell{:}, uCell{:});
   end

end