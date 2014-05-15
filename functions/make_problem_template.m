function prob = make_problem_template()

% Define global parameters here (including control bounds is a good idea)
% Include anything here that you want to be able to change at runtime
global umin umax ANY_ADDITIONAL_PARAMETERS;

% This block attaches the functions defined below as fields to the problem
% structure. If some are not needed delete them here and below.
prob.g = @g;
prob.dgdx = @dgdx;
prob.dgdu = @dgdu;

prob.f = @f;
prob.dfdx_times_vec = @dfdx_times_vec;
prob.dfdu_times_vec = @dfdu_times_vec;

prob.adjointRHS = @adjointRHS; 
prob.dHdu = @dHdu;
prob.ControlChar = @ControlChar; 
prob.optJac = @optJac;
prob.ControlBounds = [umin umax]; % lower bound in 1st column, upper bound in 2nd
                                  % additional rows for more controls if needed

%----------------------------------------------------------------------------                               
% Subfunctions
%----------------------------------------------------------------------------
% All functions except optJac should be vectorized in the second component
% for example:  t + x(1, :).*u(2, :) etc. to handle t = row-vector.

   function value = g(t, x, u)
      value = ;
   end

   function value = dgdx(t, x, u)
      value = ;
   end

   function value = dgdu(t, x, u)
      value = ;
   end

   function value = f(t, x, u)
      value = ;
   end

   function value = dfdx_times_vec(t, x, u, v)
      value = ;
   end

   function value = dfdu_times_vec(t, x, u, v)
      value = ;
   end




   function value = adjointRHS(t, x, lam, u)
      value = ;
   end

   function value = dHdu(t, x, lam, u)
      value = ;
   end

   function value = ControlChar(t, x, lam)
      value = min(umax, max(umin, CONTENT_HERE));
   end

   function J = optJac(t, y)
      nSTATES = size(y, 1)/2;
      J = zeros(2*nSTATES, 2*nSTATES);
      x = y(1:nSTATES); lam = y(nSTATES+1:end); u = ControlChar(t, x, lam);
      if (u == umin) || (u==umax)
         % Define any relevant partial derivatives of u as zero here
      else
         % Define any relevant partial derivatives of u according to the 
         % derivative of ControlChar here
      end
      
      % Give formulas for each entry of J here
   end
end