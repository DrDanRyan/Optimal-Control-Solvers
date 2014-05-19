function prob = rk4_test_problem()

% Parameters (actual values are assigned in the "solve_test_problem.m" script)
global c m r umax;

% This block attaches the functions defined below as fields to the problem
% structure.
prob.F = @F;
prob.dFdx_times_vec = @dFdx_times_vec;
prob.dFdu_times_vec = @dFdu_times_vec;
prob.ControlBounds = [0 umax]; %lower bound in 1st column, upper bound in 2nd
                               %additional rows for more controls if needed


   function value = F(t, x, u)
      value = [x.*(m - x) - u;
               exp(-r*t).*(x.^2 + c*u.^2)];
   end

   function value = dFdx_times_vec(t, x, ~, v)
      value = [(m - 2*x).*v(1,:) + 2*exp(-r*t).*x.*v(2,:); 
               0];
   end

   function value = dFdu_times_vec(t, ~, u, v)
      value = -v(1,:) + 2*c*exp(-r*t).*u.*v(2,:);
   end

end