classdef PWConstantControl < Control
   % A piecewise constant control function
   
   properties
      nControlIntervals
      nControls
   end
   
   methods
      
      function [Lb, Ub] = compute_nlp_bounds(obj, controlBounds)
         Lb = controlBounds(:,1)*ones(1, obj.nControlPts);
         Lb = reshape(Lb, [], 1);

         Ub = controlBounds(:,2)*ones(1, obj.nControlPts);
         Ub = reshape(Ub, [], 1);
      end     
      
      function dJdv = compute_dJdv(obj, dJdu)
         
      end
      
      function u = compute_u(obj, v)
         
      end
   end
   
end

