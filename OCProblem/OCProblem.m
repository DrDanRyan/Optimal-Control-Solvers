classdef OCProblem < handle
   
   properties (Abstract)
      ControlBounds
   end
   
   
   methods (Abstract)
      
      dxdt = F(obj, t, y, u)
      value = dFdx_times_vec(obj, t, y, u, v)
      value = dFdu_times_vec(obj, t, y, u, v)      
      
   end
   
end

