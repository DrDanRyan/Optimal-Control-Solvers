function [x, lam, uNew] = compute_x_lam(prob, x0, tspan, u, RelTol, AbsTol)
% Auxillary function that comes up in many of the solution algorithms

T0 = tspan(1); TF = tspan(end); lam0 = 0*x0;
nStates = length(x0);

forwardRHS = @(t, x) prob.stateRHS(t, x, u(t));
[tout, xout] = odevr7(forwardRHS, [T0, TF], x0, RelTol, AbsTol);
x = vectorInterpolant(tout', xout', 'pchip');

backwardRHS = @(t, lam) prob.adjointRHS(t, x(t), lam, u(t));
[tPts, lamPts] = odevr7(backwardRHS, [TF, T0], lam0, RelTol, AbsTol);
tPts = flip(tPts');
lamPts = flip(lamPts');
uPts = prob.ControlChar(tPts, x(tPts), lamPts);

lam = vectorInterpolant(tPts, lamPts, 'pchip');
uNew = vectorInterpolant(tPts, uPts, 'pchip');
end