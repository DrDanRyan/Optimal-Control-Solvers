classdef ControlType < handle
% This defines the interface required for a ControlType subclass
   
   methods (Abstract)
      dJdv = compute_dJdv(obj, dJdu)
      u = compute_u(obj, t, v)
      v = compute_initial_v(obj, controlBounds)
      uFunc = compute_uFunc(v)
      
      % These are optional if nlp bound or nonlinear constraints are in effect
      % [Lb, Ub] = compute_nlp_bounds(obj, controlBounds)
      % [c, ceq] = compute_nonlcon(obj, v)
      
   end

   
end

