classdef OCProblem < handle
   % Defines the abstract interface for OCProblem instances
   
   properties
      x0
      T
      nSteps 
   end
   
   properties (Dependent)
      h
      tspan
      nStates
   end
   
   methods (Abstract)
      dxdt = f(t, x, u)
      dfdx = dfdx(t, x, u)
      dfdu = dfdu(t, x, u)
      
      dJdt = g(t, x, u)          
      dgdx = dgdx(t, x, u)
      dgdu = dgdu(t, x, u)
   end
   
   methods
      
      function h = get.h(obj)
         h = obj.T/obj.nSteps;
      end
      
      
      function tspan = get.tspan(obj)
         tspan = linspace(0, obj.T, obj.nSteps+1);
      end
      
      
      function nStates = get.nStates(obj)
         nStates = length(obj.x0);
      end
           
      
      function midPts = compute_midpoints(vec)
         midPts = .5*(vec(:,1:end-1) + vec(:,2:end));
      end
      
      
      function [x, J, xK] = compute_states(obj, u)
         x = zeros(obj.nStates, obj.nSteps+1);
         x(:,1) = obj.x0;
         xK = zeros(obj.nStates, obj.nSteps, 3);
         cumObj = zeros(1, obj.nSteps+1);
         
         t = obj.tspan;
         h = obj.h;
         tHalf = obj.compute_midpoints(t);
         uHalf = obj.compute_midpoints(u);
         
         for i = 1:obj.nSteps
            y1 = obj.f(t(i), x(:,i), u(:,i));
            xK(:,i,1) = x(:,i) + .5*h*y1;
            
            y2 = obj.f(tHalf(i), xK(:,i,1), uHalf(:,i));
            xK(:,i,2) = x(:,i) + .5*h*y2;
            
            y3 = obj.f(tHalf(i), xK(:,i,2), uHalf(:,i));
            xK(:,i,3) = x(:,i) + h*y3;
            
            y4 = obj.f(t(i+1), xK(:,i,3), u(:,i+1));
            
            x(:,i+1) = x(:,i) + h*(y1 + 2*y2 + 2*y3 + y4)/6;
            
            y1 = obj.g(t(i), x(:,i), u(:,i));
            y2 = obj.g(tHalf(i), xK(:,i,1), uHalf(:,i));
            y3 = obj.g(tHalf(i), xK(:,i,2), uHalf(:,i));
            y4 = obj.g(t(i+1), xK(:,i,3), u(:,i+1));
            
            cumObj(i+1) = cumObj(i) + h*(y1 + 2*y2 + 2*y3 + y4)/6;
         end
         J = cumObj(end);         
      end
      
      function lam = compute_adjoints(obj, x, xK, u)
         lam = zeros(obj.nStates, obj.nSteps+1);
         t = obj.tspan;
         h = obj.h;
         tHalf = obj.compute_midpoints(t);
         uHalf = obj.compute_midpoints(u);
         
         for i = obj.nSteps:-1:1
            y1 = obj.dfdx(t(i), x(:,i), u(:,i));
            y2 = obj.dfdx(tHalf(i), xK(:,i,1), uHalf(:,i))*...
               (eye(obj.nStates) + .5*h*y1);  % TODO check if matrix order is correct in multiplication
            y3 = obj.dfdx(tHalf(i), xK(:,i,2), uHalf(:,i))*...
               (eye(obj.nStates) + .5*h*y2);
            y4 = obj.dfdx(t(i+1), xK(:,i,3), u(:,i+1))*...
               (eye(obj.nStates) + h*y3);
            
            
         end
      end
      
      function dJdu = compute_objective_gradient(x, lam, u)
         
      end
   end
      
end

