function x = compute_states(prob, x0, t, u)

nSTATES = length(x0);
nTSTEPS = length(t);
h = t(2) - t(1); % t should be equally spaced TODO: raise exception if not true

x = zeros(nSTATES, nTSTEPS);
x(:,1) = x0;

tHalf = t(1:end-1) + h/2;
uHalf = .5*(u(:,1:end-1) + u(:,2:end));

for i = 1:nTSTEPS-1
   k1 = prob.stateRHS(t(i), x(:,i), u(:,i));
   k2 = prob.stateRHS(tHalf(i), x(:,i) + .5*h*k1, uHalf(:,i));
   k3 = prob.stateRHS(tHalf(i), x(:,i) + .5*h*k2, uHalf(:,i));
   k4 = prob.stateRHS(t(i+1), x(:,i) + h*k3, uHalf(:,i+1));
   x(:,i+1) = x(:,i) + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end


end

