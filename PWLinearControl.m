classdef PWLinearControl < ControlType

   properties
      controlBounds % [Lb, Ub]
      uspan
   end
   
   
   methods
      function obj = PWLinearControl(nCONTROL_PTS, timeRange, controlBounds)
         obj.uspan = linspace(timeRange(1), timeRange(2), nCONTROL_PTS);
         Lb = controlBounds(:,1)*ones(1, nCONTROL_PTS);
         Lb = reshape(Lb, [], 1);
         
         Ub = controlBounds(:,2)*ones(1, nCONTROL_PTS);
         Ub = reshape(Ub, [], 1);
         obj.controlBounds = [Lb, Ub];
      end
      
      function dJdv = compute_dJdv(obj, dJdu)
      
      end
      
      
      function u = compute_u(obj, t, v)
         uFunc = obj.compute_uFunc(v);
         u = uFunc(t);
      end
      
      
      function v = compute_initial_v(obj, vSize)
         if isempty(obj.controlBounds)
            v = zeros(vSize, 1);
         else
            v = max(obj.controlBounds(:,1), 0);
         end
      end
      
      
      function [c, ceq] = compute_nonlcon(~, ~)
         c = [];
         ceq = [];
      end
      
      
      function uFunc = compute_uFunc(obj, v)
         uFunc  = vectorInterpolant(obj.uspan, v, 'linear');
      end
   end
   
end

