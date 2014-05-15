classdef RK4Integrator < Integrator
   
   properties
      t   % tspan augmented with midpoints of each interval
      h   % diff(tspan)
      nSTEPS   % number of integration steps = length(h)
      
      % Stores state values and intermediate state estimates computed during 
      % the last forward pass. These values are used to compute the adjoints.
      xK
   end
   
   
   methods
      function obj = RK4Integrator(tspan)
         obj.h = diff(tspan);
         obj.nSTEPS = length(obj.h);
         
         % t is a time vector that also includes the interval midpoints
         t = zeros(1, 2*obj.nSTEPS+1); 
         t(1:2:end) = tspan;
         t(2:2:end-1) = (tspan(1:end-1) + tspan(2:end))/2;
         obj.t = t;
      end
      
      
      function [x, J] = compute_states(obj, prob, x0, u)
         nSTATES = length(x0) + 1;
         
         % store x, xk1, xk2, and xk3 at each time pt
         obj.xK = nan([nSTATES, obj.nSTEPS+1, 4]); 
         obj.xK(:,1,1) = [x0; 0];

         for i = 1:obj.nSTEPS
            % Perform single RK-4 Step
            % F1, F2, F3, F4 refer to the RK approx values of the state rhs

            F1 = prob.F(obj.t(2*i-1), obj.xK(:,i,1), u(:,2*i-1));
            obj.xK(:,i,2) = obj.xK(:,i,1) + obj.h(i)/2*F1;

            F2 = prob.F(obj.t(2*i), obj.xK(:,i,2), u(:,2*i));
            obj.xK(:,i,3) = obj.xK(:,i,1) + obj.h(i)/2*F2;

            F3 = prob.F(obj.t(2*i), obj.xK(:,i,3), u(:,2*i));
            obj.xK(:,i,4) = obj.xK(:,i,1) + obj.h(i)*F3;

            F4 = prob.F(obj.t(2*i+1), obj.xK(:,i,4), u(:,2*i+1));

            obj.xK(:,i+1,1) = obj.xK(:,i,1) ...
                                 + obj.h(i)/6*(F1 + 2*F2 + 2*F3 + F4);         
         end         

         x = obj.xK(:,:,1);
         J = x(end,end);
      end
      
      
      function [lam, dJdu] = compute_adjoints(obj, prob, u)
         
         nSTATES = size(obj.xK, 1);
         lam = zeros(nSTATES, obj.nSTEPS+1);      
         lam(end,end) = 1;
         dJdk = nan(nSTATES, obj.nSTEPS, 4);

         for i = obj.nSTEPS:-1:1
            dJdk(:,i,4) = obj.h(i)/6*lam(:,i+1);
            dJdx3 = prob.dFdx_times_vec(obj.t(2*i+1), obj.xK(:,i,4), ...
                                        u(:,2*i+1), dJdk(:,i,4));

            dJdk(:,i,3) = obj.h(i)/3*lam(:,i+1) + obj.h(i)*dJdx3;
            dJdx2 = prob.dFdx_times_vec(obj.t(2*i), obj.xK(:,i,3), u(:,2*i), ...
                                        dJdk(:,i,3));

            dJdk(:,i,2) = obj.h(i)/3*lam(:,i+1) + obj.h(i)/2*dJdx2;
            dJdx1 = prob.dFdx_times_vec(obj.t(2*i), obj.xK(:,i,2), u(:,2*i), ...
                                        dJdk(:,i,2));

            dJdk(:,i,1) = obj.h(i)/6*lam(:,i+1) + obj.h(i)/2*dJdx1;
            lam(:,i) = lam(:,i+1) + dJdx1 + dJdx2 + dJdx3 + ...
                  prob.dFdx_times_vec(obj.t(2*i-1), obj.xK(:,i,1), ...
                                      u(:,2*i-1), dJdk(:,i,1));
         end  

         if nargout > 1
            dJdu = obj.compute_dJdu(prob, u, dJdk);
         end
      end
      
      
      function dJdu = compute_dJdu(obj, prob, u, dJdk)
         dJdu = zeros(size(u)); % ~ nCONTROLS x 2*nSTEPS+1
         
         % Left end point
         dJdu(:,1) = prob.dFdu_times_vec(obj.t(1), obj.xK(:,1,1), u(:,1), ...
                                         dJdk(:,1,1));
         
         % RK step interval mid points
         dJdu(:,2:2:end-1) = ...
            prob.dFdu_times_vec(obj.t(2:2:end-1), obj.xK(:,1:end-1,2), ...
                                u(:,2:2:end-1), dJdk(:,:,2)) ...
            + prob.dFdu_times_vec(obj.t(2:2:end-1), obj.xK(:,1:end-1,3), ...
                                  u(:,2:2:end-1), dJdk(:,:,3));
         
         % RK step interval endpts (except for left and right of full interval
         dJdu(:,3:2:end-2) = ...
            prob.dFdu_times_vec(obj.t(3:2:end-2), obj.xK(:,2:end-1,1), ...
                                u(:,3:2:end-2), dJdk(:,2:end,1)) ...
            + prob.dFdu_times_vec(obj.t(3:2:end-2), obj.xK(:,1:end-2,4), ...
                                  u(:,3:2:end-2), dJdk(:,1:end-1,4));
         
         % Right end point
         dJdu(:,end) = prob.dFdu_times_vec(obj.t(end), obj.xK(:,end-1,4), ...
                                           u(:,end), dJdk(:,end,4)); 
      end
      
   end
end

