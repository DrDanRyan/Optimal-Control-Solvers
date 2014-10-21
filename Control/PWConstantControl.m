classdef PWConstantControl < Control
   % A piecewise constant control function
   
   properties
      nControlIntervals
      nControls
      intervalStarts
      B % basis matrix
   end
   
   methods
      function obj = PWConstantControl(t, nControlIntervals, nControls)
         obj.nControlIntervals = nControlIntervals;
         obj.nControls = nControls;
         obj.intervalStarts = linspace(t(1), t(end), nControlIntervals+1);
         obj.intervalStarts = obj.intervalStarts(1:end-1);
         obj.B = obj.compute_basis_matrix(t);
      end
      
      
      function [Lb, Ub] = compute_nlp_bounds(obj, controlBounds)
         Lb = controlBounds(:,1)*ones(1, obj.nControlIntervals);
         Lb = reshape(Lb, [], 1);

         Ub = controlBounds(:,2)*ones(1, obj.nControlIntervals);
         Ub = reshape(Ub, [], 1);
      end     
      
      
      function B = compute_basis_matrix(obj, t)
         B = zeros(obj.nControlIntervals, length(t));
         for i=1:obj.nControlIntervals-1
            B(i,:) = double(t >= obj.intervalStarts(i) & t < obj.intervalStarts(i+1));
         end
         
         B(end,:) = double(t >= obj.intervalStarts(end));
      end
      
      
      function dJdv = compute_dJdv(obj, dJdu)
         dJdv = dJdu*obj.B';
         dJdv = reshape(dJdv, [], 1);
      end
      
      
      function u = compute_u(obj, v)
         v = reshape(v, obj.nControls, []);
         u = v*obj.B;
      end
      
      
      function v = compute_initial_v(obj, u0)
         v = repmat(u0, obj.nControlIntervals, 1);
      end
      
      
      function uFunc = compute_uFunc(obj, v)
         v = reshape(v, obj.nControls, []);
         uFunc  = vectorInterpolant(obj.intervalStarts, v, 'previous');
      end
      
   end
   
end

