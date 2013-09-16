function [x, lam, J] = compute_x_lam_J(prob, x0, tspan, u, RelTol, AbsTol)
% Auxillary function that comes up in solve_bvp

T0 = tspan(1); TF = tspan(end); lam0 = 0*x0; nStates = length(x0);

   function value = forwardRHS(t, y)
      tLENGTH = length(t);
      value = zeros(nStates+1, tLENGTH);
      xVec = y(1:nStates, :);
      value(1:nStates, :) = prob.stateRHS(t, xVec, u(t));
      value(end, :) = prob.objective(t, xVec, u(t));
   end

[tout, xAugout] = odevr7(@forwardRHS, [T0, TF], [x0; 0], RelTol, AbsTol);
J = xAugout(end, end);
x = vectorInterpolant(tout', xAugout(:, 1:nStates)', 'pchip');

backwardRHS = @(t, lam) prob.adjointRHS(t, x(t), lam, u(t));
[tout, lamout] = odevr7(backwardRHS, [TF, T0], lam0, RelTol, AbsTol);
lam = vectorInterpolant(flip(tout', 2), flip(lamout', 2), 'pchip');
end