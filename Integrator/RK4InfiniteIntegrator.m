classdef RK4InfiniteIntegrator < Integrator
   
   properties
      t
      integrator1
      integrator2
      uStar
   end
   
   
   methods
      function obj = RK4InfiniteIntegrator(tspan, tspanExtra, uStar)
         obj.integrator1 = RK4Integrator(tspan);
         obj.integrator2 = RK4Integrator(tspanExtra);
         obj.uStar = uStar*ones(size(obj.integrator2.t)); % col vec times row vec
         obj.t = obj.integrator1.t; % recode as Dependent property?
      end
      
      
      function [x, J] = compute_states(obj, prob, x0, u)
         [x, J1] = obj.integrator1.compute_states(prob, x0, u);
         [~, J2] = obj.integrator2.compute_states(prob, x(1:end-1,end), obj.uStar);
         J = J1 + J2;
      end
      
      
      function [lam, dJdu] = compute_adjoints(obj, prob, u)
         lam2 = obj.integrator2.compute_adjoints(prob, obj.uStar);
         [lam, dJdu] = obj.integrator1.compute_adjoints(prob, u, lam2(:,1));
      end
   end
   
end

