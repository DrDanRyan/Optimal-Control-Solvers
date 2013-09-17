function prob = make_test_problem()

% Parameters (actual values are assigned in the "solve_test_problem.m" script)
global c m umax;

% This block attaches the functions defined below as fields to the problem
% structure.
prob.objective = @objective;
prob.stateRHS = @stateRHS;
prob.adjointRHS = @adjointRHS; % dLam/dt = -dH/dx
prob.dHdu = @dHdu;
prob.ControlChar = @ControlChar; %Found by applying Pontryagin
prob.optJac = @optJac;
prob.ControlBounds = [0 umax]; %lower bound in 1st column, upper bound in 2nd
                               %additional rows for more controls if needed

   function value = objective(t, x, u)
      value = x.^2 + c*u.^2;
   end

   function value = stateRHS(t, x, u)
      value = x.*(m - x) - u;
   end

   function value = adjointRHS(t, x, lam, u)
      value = -2*x - m*lam + 2*lam.*x;
   end

   function value = dHdu(t, x, lam, u)
      value = 2*c*u - lam;
   end

   function value = ControlChar(t, x, lam)
      value = min(umax, max(0, lam/(2*c)));
   end

   function J = optJac(t, y)
      J = zeros(2, 2);
      x = y(1); lam = y(2); u = ControlChar(t, x, lam);
      if (u == 0) || (u==umax)
         dudLam = 0;
      else
         dudLam = 1/(2*c);
      end
      J(1, 1) = m - 2*x;
      J(1, 2) = -dudLam;
      J(2, 1) = -2 + 2*lam;
      J(2, 2) = -m + 2*x;
   end

end