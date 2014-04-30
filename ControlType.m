classdef ControlType < handle
% This defines the interface required for a ControlType subclass

   properties
      controlBounds
      isConstraintDerivative = false
   end
   
   methods (Abstract)
      dJdv = compute_dJdv(obj, dJdu)
      u = compute_u(obj, t, v)
      v = compute_initial_v(obj)
      [c, ceq] = compute_nonlcon(obj, v)
      uFunc = compute_uFunc(v)
   end

   
end

