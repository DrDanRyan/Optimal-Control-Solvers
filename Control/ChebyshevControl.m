classdef ChebyshevControl < Control

   properties
      nControlBasis
      nControls     
      controlPts
      
      B
   end
   
   
   methods
      function obj = ChebyshevControl(t, nControlBasis, nControls)
         obj.nControlBasis = nControlBasis;
         obj.nControls = nControls;
         obj.controlPts = linspace(t(1), t(end), nControlBasis);
         obj.B = obj.compute_basis_matrix(t);
      end
      
      
      function B = compute_basis_matrix(obj, t)
         B = zeros(obj.nControlBasis, length(t));
         tTrans = 2*(t - t(1))/(t(end) - t(1)) - 1; % tTrans in [-1, 1]
         
         B(1,:) = 1;
         B(2,:) = tTrans;
         
         for i = 3:obj.nControlBasis
            B(i,:) = 2*tTrans.*B(i-1,:) - B(i-2,:);
         end
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
         v = [u0; zeros((obj.nControlBasis-1)*obj.nControls, 1)];
      end
      
      
      function [A, b] = compute_lincon(obj, ControlBounds)
         
      end


   end
   
end

