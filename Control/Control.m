classdef Control < handle
% This defines the interface for a Control object
   
   methods (Abstract)
      dJdv = compute_dJdv(obj, dJdu)
      u = compute_u(obj, t, v)
      v = compute_initial_v(obj, u0)
      uFunc = compute_uFunc(v)
      
      % These are optional if nlp bound or nonlinear constraints are in effect
      % [Lb, Ub] = compute_nlp_bounds(obj, controlBounds)
      % [c, ceq] = compute_nonlcon(obj, v)
      
   end

   
end

