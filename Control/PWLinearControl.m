classdef PWLinearControl < Control

   properties
      nControlPts
      nControls     
      controlPts
      
      B  % basis functions evaluated at points required by integrator (t)
   end
   
   
   methods
      function obj = PWLinearControl(t, nControlPts, nControls)
         obj.nControlPts = nControlPts;
         obj.nControls = nControls;
         obj.controlPts = linspace(t(1), t(end), nControlPts);
         obj.B = obj.compute_basis_matrix(t);
      end
      
      
      function [Lb, Ub] = compute_nlp_bounds(obj, controlBounds)
         
         Lb = controlBounds(:,1)*ones(1, obj.nControlPts);
         Lb = reshape(Lb, [], 1);
         
         Ub = controlBounds(:,2)*ones(1, obj.nControlPts);
         Ub = reshape(Ub, [], 1);
      end
      
      
      function B = compute_basis_matrix(obj, t)
         B = zeros(obj.nControlPts, length(t));
         
         % Basis function at left end of time interval
         tent = griddedInterpolant(obj.controlPts(1:2), [1, 0], ...
                                   'linear', 'nearest');
         B(1,:) = tent(t);
         
         % All basis functions in the interior
         for i = 2:obj.nControlPts-1
            tent = griddedInterpolant(obj.controlPts(i-1:i+1), [0, 1, 0], ...
                                      'linear', 'nearest');
            B(i,:) = tent(t);            
         end
         
         % Basis function at right end of time interval          
         tent = griddedInterpolant(obj.controlPts(end-1:end), [0, 1], ...
                                   'linear', 'nearest');
         B(end,:) = tent(t);
      end
      
      
      function dJdv = compute_dJdv(obj, dJdu)
         dJdv = dJdu*obj.B';
         dJdv = reshape(dJdv, [], 1);
      end
      
      
      function u = compute_u(obj, v)
         v = reshape(v, obj.nControls, []);
         u = v*obj.B;
      end
      
      
      function v = compute_initial_v(obj, boundMatrix)
         % boundMatrix ~ nControls x 2
         v = max(boundMatrix(:,1), 0);
         v = repmat(v, obj.nControlPts, 1);
      end
      
      
      function uFunc = compute_uFunc(obj, v)
         v = reshape(v, obj.nControls, []);
         uFunc  = vectorInterpolant(obj.controlPts, v, 'linear');
      end

   end
   
end

