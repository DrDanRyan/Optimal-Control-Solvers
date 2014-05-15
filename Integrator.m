classdef Integrator < handle
   % This defines the Integrator interface. An Integrator's main responsibility
   % is to compute the state, adjoint, functional and functional gradient given
   % control values at a specified time grid.
   
   properties (Abstract)
      % a vector of time pts where the control value has an effect 
      % on the objective
      t 
   end
   
   methods (Abstract)
      [x, J] = compute_states(obj, prob, x0, u)
      [lam, dJdu] = compute_adjoints(obj, prob, u)
   end
   
end

