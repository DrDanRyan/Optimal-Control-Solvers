classdef TestOCProblem < OCProblem
   
   properties
      ControlBounds
      
      c
      m
      r
   end
   
   
   methods
      function obj = TestOCProblem(p, ControlBounds)
         obj.ControlBounds = ControlBounds;
         
         obj.c = p.c;
         obj.m = p.m;
         obj.r = p.r;
      end
      
      
      function value = F(obj, t, y, u)
         x = y(1,:);
         value = [x.*(obj.m - x) - u;
                  exp(-obj.r*t).*(x.^2 + obj.c*u.^2)];
      end
      
      
      function value = dFdx_times_vec(obj, t, y, ~, v)
         x = y(1,:);
         value = (obj.m - 2*x).*v(1,:) + 2*exp(-obj.r*t).*x.*v(2,:);
      end
      
      
      function value = dFdu_times_vec(obj, t, ~, u, v)
         value = -v(1,:) + 2*obj.c*exp(-obj.r*t).*u.*v(2,:);
      end
      
   end
end

