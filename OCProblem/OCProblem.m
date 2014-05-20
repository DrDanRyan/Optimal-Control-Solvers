classdef OCProblem < handle
   
   properties (Abstract)
      ControlBounds
   end
   
   
   methods (Abstract)
      
      % This function computes the state equations RHS. The objective integrand
      % should be the formula in the last row.
      dxdt = F(obj, t, y, u)
      
      % This function computes dFdx'*v. Because the objective variable, y(end),
      % is not used to compute F, the last row of value will be zero.
      value = dFdx_times_vec(obj, t, y, u, v)
      
      % This function computes dFdu'*v.
      value = dFdu_times_vec(obj, t, y, u, v)      
      
   end
   
end

