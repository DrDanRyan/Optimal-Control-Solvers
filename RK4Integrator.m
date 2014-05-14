classdef RK4Integrator < Integrator
   
   properties
      t
   end
   
   
   methods
      function obj = RK4Integrator(t)
         obj.t = t;
      end
      
      
      function [x, J] = compute_states(obj, prob, x0, isSave)
         
      end
      
      
      function [lam, dJdu] = compute_adjoints(obj, prob, x, u)
         
      end
      
   end
end

