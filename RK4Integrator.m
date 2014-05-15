classdef RK4Integrator < Integrator
   
   properties
      t   % tspan augmented with midpoints of each interval
   end
   
   
   methods
      function obj = RK4Integrator(tspan)
         nSTEPS = length(tspan - 1);
         t = zeros(1, 2*nSTEPS+1); % a time vector that also includes the interval midpoints
         t(1:2:end) = tspan;
         t(2:2:end-1) = (tspan(1:end-1) + tspan(2:end))/2;
         obj.t = t;
      end
      
      
      function [x, J] = compute_states(obj, prob, x0, u)
%          x = nan([STATE_SHAPE, 4]); % store x, xk1, xk2, and xk3 at each time pt
%          x(:,1,1) = [x0; 0];
% 
%          for i = 1:nSTEPS
%             % Perform single RK-4 Step
%             % F1, F2, F3, F4 refer to the RK approx values of the state rhs
% 
%             F1 = prob.F(t(2*i-1), x(:,i,1), u(:,2*i-1));
%             x(:,i,2) = x(:,i,1) + h(i)/2*F1;
% 
%             F2 = prob.F(t(2*i), x(:,i,2), u(:,2*i));
%             x(:,i,3) = x(:,i,1) + h(i)/2*F2;
% 
%             F3 = prob.F(t(2*i), x(:,i,3), u(:,2*i));
%             x(:,i,4) = x(:,i,1) + h(i)*F3;
% 
%             F4 = prob.F(t(2*i+1), x(:,i,4), u(:,2*i+1));
% 
%             x(:,i+1,1) = x(:,i,1) + h(i)/6*(F1 + 2*F2 + 2*F3 + F4);         
%          end         
% 
%          J = x(end,end,1);
      end
      
      
      function [lam, dJdu] = compute_adjoints(obj, prob, x, u)
%          lam = zeros(STATE_SHAPE);      
%          lam(end,end) = 1;
% 
%          dJdk = nan(nSTATES, nSTEPS, 4);
% 
%          for i = nSTEPS:-1:1
%             dJdk(:,i,4) = h(i)/6*lam(:,i+1);
%             dJdx3 = prob.dFdx_times_vec(t(2*i+1), x(:,i,4), u(:,2*i+1), ...
%                                         dJdk(:,i,4));
% 
%             dJdk(:,i,3) = h(i)/3*lam(:,i+1) + h(i)*dJdx3;
%             dJdx2 = prob.dFdx_times_vec(t(2*i), x(:,i,3), u(:,2*i), ...
%                                         dJdk(:,i,3));
% 
%             dJdk(:,i,2) = h(i)/3*lam(:,i+1) + h(i)/2*dJdx2;
%             dJdx1 = prob.dFdx_times_vec(t(2*i), x(:,i,2), u(:,2*i), ...
%                                         dJdk(:,i,2));
% 
%             dJdk(:,i,1) = h(i)/6*lam(:,i+1) + h(i)/2*dJdx1;
%             lam(:,i) = lam(:,i+1) + dJdx1 + dJdx2 + dJdx3 + ...
%                   prob.dFdx_times_vec(t(2*i-1), x(:,i,1), u(:,2*i-1), dJdk(:,i,1));
%          end  
% 
%          if nargout > 1
%             dJdu = compute_dJdu(x, u, dJdk);
%             dJdv = control.compute_dJdv(dJdu);
%          end
      end
      
      function dJdu = compute_dJdu(obj, x, u, dJdk)
%          dJdu = zeros(size(u)); % ~ nCONTROLS x 2*nSTEPS+1
%          
%          % Left end point
%          dJdu(:,1) = prob.dFdu_times_vec(t(1), x(:,1,1), u(:,1), dJdk(:,1,1));
%          
%          % RK step interval mid points
%          dJdu(:,2:2:end-1) = ...
%             prob.dFdu_times_vec(t(2:2:end-1), x(:,1:end-1,2), ...
%                                 u(:,2:2:end-1), dJdk(:,:,2)) ...
%             + prob.dFdu_times_vec(t(2:2:end-1), x(:,1:end-1,3), ...
%                                   u(:,2:2:end-1), dJdk(:,:,3));
%          
%          % RK step interval endpts (except for left and right of full interval
%          dJdu(:,3:2:end-2) = ...
%             prob.dFdu_times_vec(t(3:2:end-2), x(:,2:end-1,1), ...
%                                 u(:,3:2:end-2), dJdk(:,2:end,1)) ...
%             + prob.dFdu_times_vec(t(3:2:end-2), x(:,1:end-2,4), ...
%                                   u(:,3:2:end-2), dJdk(:,1:end-1,4));
%          
%          % Right end point
%          dJdu(:,end) = prob.dFdu_times_vec(t(end), x(:,end-1,4), ...
%                                            u(:,end), dJdk(:,end,4)); 
      end
      
   end
end

