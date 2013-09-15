function J = compute_J(prob, x0, tspan, u, RelTol, AbsTol)
% Auxillary function that comes up in fb_sweep

T0 = tspan(1); TF = tspan(end); nSTATES = length(x0);

   function value = forwardRHS(t, y)
      value = zeros(nSTATES+1, size(t, 2));
      xVec = y(1:nSTATES, :);
      uVec = u(t);
      value(1:nSTATES, :) = prob.stateRHS(t, xVec, uVec);
      value(end, :) = prob.objective(t, xVec, uVec);
   end

[~, xAugout] = odevr7(@forwardRHS, [T0, TF], [x0; 0], RelTol, AbsTol);
J = xAugout(end, end);
end