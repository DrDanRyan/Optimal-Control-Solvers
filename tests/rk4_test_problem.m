function prob = rk4_test_problem()

% Parameters (actual values are assigned in the "solve_test_problem.m" script)
global c m umax;

% This block attaches the functions defined below as fields to the problem
% structure.
prob.g = @g;
prob.f = @f;
prob.dgdx = @dgdx;
prob.dgdu = @dgdu;
prob.dfdx_times_vec = @dfdx_times_vec;
prob.dfdu_times_vec = @dfdu_times_vec;
prob.ControlBounds = [0 umax]; %lower bound in 1st column, upper bound in 2nd
                               %additional rows for more controls if needed

   function value = g(t, x, u)
      value = x.^2 + c*u.^2;
   end

   function value = dgdx(t, x, u)
      value = 2*x;
   end

   function value = dgdu(t, x, u)
      value = 2*c*u;
   end

   function value = f(t, x, u)
      value = x.*(m - x) - u;
   end

   function value = dfdx_times_vec(t, x, u, v)
      value = (m - 2*x).*v;
   end

   function value = dfdu_times_vec(t, x, u, v)
      value = -v;
   end

end