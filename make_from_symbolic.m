function prob = make_from_symbolic(symObjective, symStateRHS, nStates, nControls, ...
                                    params, bounds)
param_names = fields(params);
param_vals = struct2cell(params);
t = sym('t');
x = sym('x', [1 nStates]);
lam = sym('lam', [1 nStates]);
u = sym('u', [1 nControls]);

% Hamiltonian
symH(t, x, lam, u) = symObjective + lam*symStateRHS;
 
% % Adjoint ODE rhs
symAdjointRHS = -gradient(symH, x);

% dH/du
symdHdu = gradient(symH, u);

% Optimality condition
uTemp = unpack(u);
ControlCharTemp = cell(1, nControls);
[ControlCharTemp{:}] = solve(symdHdu, uTemp{:});
symControlChar(t, x, lam) = transpose(sym(ControlCharTemp));
clear uTemp ControlCharTemp

% Convert symbolic function to matlab function handles
objectiveTemp = matlabFunction(symObjective, 'vars', [t, x, u, param_names']);
stateRHSTemp = matlabFunction(symStateRHS, 'vars', [t, x, u, param_names']);
adjointRHSTemp = matlabFunction(symAdjointRHS, 'vars', [t, x, lam, u, param_names']);
dHduTemp = matlabFunction(symdHdu, 'vars', [t, x, lam, u, param_names']);
ControlCharTemp = matlabFunction(symControlChar, 'vars', [t, x, lam, param_names']);

prob.objective = @objective;
prob.stateRHS = @stateRHS;
prob.adjointRHS = @adjointRHS;
prob.dHdu = @dHdu;
prob.ControlChar = @ControlChar;
prob.ControlBounds = bounds;
umin = bounds(:,1);
umax = bounds(:,2);

   function value = objective(t, x, u)
      xCell = cell(1, nStates);
      uCell = cell(1, nControls);
      for i = 1:nStates
         xCell{i} = x(i,:);
      end

      for i = 1:nControls
         uCell{i} = u(i,:);
      end

      value = objectiveTemp(t, xCell{:}, uCell{:}, param_vals{:});
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

      value = stateRHSTemp(t, xCell{:}, uCell{:}, param_vals{:});
   end

   function value = adjointRHS(t, x, lam, u)
      xCell = cell(1, nStates);
      lamCell = cell(1, nStates);
      uCell = cell(1, nControls);
      for i = 1:nStates
         xCell{i} = x(i, :);
         lamCell{i} = lam(i,:);
      end

      for i = 1:nControls
         uCell{i} = u(i,:);
      end

      value = adjointRHSTemp(t, xCell{:}, lamCell{:}, uCell{:}, param_vals{:});
   end

   function value = dHdu(t, x, lam, u)
      xCell = cell(1, nStates);
      lamCell = cell(1, nStates);
      uCell = cell(1, nControls);
      for i = 1:nStates
         xCell{i} = x(i, :);
         lamCell{i} = lam(i,:);
      end

      for i = 1:nControls
         uCell{i} = u(i,:);
      end

      value = dHduTemp(t, xCell{:}, lamCell{:}, uCell{:}, param_vals{:});
   end

   function value = ControlChar(t, x, lam)
      xCell = cell(1, nStates);
      lamCell = cell(1, nStates);
      for i = 1:nStates
         xCell{i} = x(i,:);
         lamCell{i} = lam(i,:);
      end

      value = ControlCharTemp(t, xCell{:}, lamCell{:}, param_vals{:});
      value = min(umax, max(umin, value)); 
   end
end

