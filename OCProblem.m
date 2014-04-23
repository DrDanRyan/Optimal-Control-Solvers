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
      
      
      function [x, J] = compute_states(obj, u)
         x = zeros(obj.nStates, obj.nSteps+1);
         x(:,1) = obj.x0;
         cumObj = zeros(1, obj.nSteps+1);
         
         t = obj.tspan;
         h = obj.h;
         tHalf = obj.compute_midpoints(t);
         uHalf = obj.compute_midpoints(u);
         
         for i = 1:obj.nSteps
            xk1 = obj.f(t(i), x(:,i), u(:,i));
            xk2 = obj.f(tHalf(i), x(:,i) + .5*h*xk1, uHalf(:,i));
            xk3 = obj.f(tHalf(i), x(:,i) + .5*h*xk2, uHalf(:,i));
            xk4 = obj.f(t(i+1), x(:,i) + h*xk3, u(:,i+1));
            x(:,i+1) = x(:,i) + h*(xk1 + 2*xk2 + 2*xk3 + xk4)/6;
            
            objk1 = obj.g(t(i), x(:,i), u(:,i));
            objk2 = obj.g(tHalf(i), x(:,i) + .5*h*xk1, uHalf(:,i));
            objk3 = obj.g(tHalf(i), x(:,i) + .5*h*xk2, uHalf(:,i));
            objk4 = obj.g(t(i+1), x(:,i) + h*xk3, u(:,i+1));
            cumObj(i+1) = cumObj(i) + h*(objk1 + 2*objk2 + 2*objk3 + objk4)/6;
         end
         J = cumObj(end);         
      end
      
      function lam = compute_adjoints(obj, u)
         
      end
   end
      
end

