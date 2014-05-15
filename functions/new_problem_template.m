function prob = new_problem_template()

% Define global parameters here (including control bounds is a good idea)
% Include anything here that you want to be able to change at runtime
global umin umax ANY_ADDITIONAL_PARAMETERS;

% This block attaches the functions defined below as fields to the problem
% structure. If some are not needed delete them here and below.
prob.F = @F;
prob.dFdx_times_vec = @dFdx_times_vec;
prob.dFdu_times_vec = @dFdu_times_vec;

prob.ControlBounds = [umin umax]; % lower bound in 1st column, upper bound in 2nd
                                  % additional rows for more controls if needed

%----------------------------------------------------------------------------                               
% Subfunctions
%----------------------------------------------------------------------------
% All functions except optJac should be vectorized in the second component
% for example:  t + x(1, :).*u(2, :) etc. to handle t = row-vector.

   function value = F(t, x, u)
      % This is the RHS of the augmented state equation (the objective value has
      % been concatenated as the last state; hence the objective integrand
      % should be the last row of F), and x(end,:) should not appear in any
      % formula in F.
      
   end

   function value = dFdx_times_vec(t, x, u, v)
      % Computes (dF/dx)'*v. Because J = x(end,:) does not appear in F, the last
      % row returned by this function should be all zeros.
      % 
      % (dF/dx)'*v = [dF/dx1'*v;
      %               dF/dx2'*v;
      %               ...
      %               dF/dxN'*v;
      %               0 ];
   end

   function value = dFdu_times_vec(t, x, u, v)
      % Computes (dF/du)'*v. 
      %
      % If u has two components:
      % (dF/du)'*v = [dF/du1'*v;
      %               dF/du2'*v];
   end
end