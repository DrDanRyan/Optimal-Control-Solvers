function prob = make_problem_template()

% Define the number of state components (used to compute optRHS in terms of
% stateRHS and adjointRHS.
nSTATES = ;

% Define global parameters here (including control bounds is a good idea)
% Include anything here that you want to be able to change at runtime
global umin umax ANY_ADDITIONAL_PARAMETERS;

% This block attaches the functions defined below as fields to the problem
% structure. If some are not needed delete them here and below.
prob.objective = @objective;
prob.stateRHS = @stateRHS;
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

   function value = objective(t, x, u)
      value = ;
   end

   function value = stateRHS(t, x, u)
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
      J = zeros(2*nSTATES, 2*nSTATES);
      x = y(1:nSTATES); lam = y(nSTATES+1:end); u = ControlChar(t, x, lam);
      if (u == 0) || (u==umax)
         % Define any relevant partial derivatives of u as zero here
      else
         % Define any relevant partial derivatives of u according to the 
         % derivative of ControlChar here
      end
      
      % Give formulas for each entry of J here
   end
end