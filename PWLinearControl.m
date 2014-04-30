classdef PWLinearControl < ControlType

   properties
      tspan
      nControlPts
      nControls     
      uspan
      dudv
   end
   
   
   methods
      function obj = PWLinearControl(tspan, nControlPts, nControls)
         obj.tspan = tspan;
         obj.nControlPts = nControlPts;
         obj.nControls = nControls;
         obj.uspan = linspace(tspan(1), tspan(end), nControlPts);       
         obj.dudv = obj.make_dudv();
      end
      
      
      function [Lb, Ub] = compute_nlp_bounds(obj, controlBounds)
         % boundMatrix ~ nCONTROLS x 2
         
         Lb = controlBounds(:,1)*ones(1, obj.nControlPts);
         Lb = reshape(Lb, [], 1);
         
         Ub = controlBounds(:,2)*ones(1, obj.nControlPts);
         Ub = reshape(Ub, [], 1);
      end
      
      
      function dudv = make_dudv(obj)
         dudv = zeros(length(obj.tspan), length(obj.uspan));
         
         tent = griddedInterpolant([obj.uspan(1), obj.uspan(2)], [1, 0], ...
                                   'linear', 'nearest');
         dudv(:,1) = tent(obj.tspan)';
         
         for i = 2:length(obj.uspan)-1
            tent = griddedInterpolant(obj.uspan(i-1:i+1), [0, 1, 0], ...
                                      'linear', 'nearest');
            dudv(:,i) = tent(obj.tspan)';
         end
         
         tent = griddedInterpolant(obj.uspan(end-1:end), [0, 1], ...
                                   'linear', 'nearest');
         dudv(:,end) = tent(obj.tspan)';
      end
      
      
      function dJdv = compute_dJdv(obj, dJdu)
         dJdv = dJdu*obj.dudv;
      end
      
      
      function u = compute_u(obj, t, v)
         uFunc = obj.compute_uFunc(v);
         u = uFunc(t);
      end
      
      
      function v = compute_initial_v(obj, boundMatrix)
         % boundMatrix ~ nControls x 2
         v = max(boundMatrix(:,1), 0);
         v = repmat(v, obj.nControlPts, 1);
      end
      
      
      function uFunc = compute_uFunc(obj, v)
         uFunc  = vectorInterpolant(obj.uspan, ...
                        reshape(v, obj.nControls, obj.nControlPts), 'linear');
      end

   end
   
end

