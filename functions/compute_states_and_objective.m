function [x, J] = compute_states_and_objective(prob, x0, t, u)

nSTATES = length(x0);
nTSTEPS = length(t);
h = t(2) - t(1); % t should be equally spaced TODO: raise exception if not true

xAug = zeros(nSTATES+1, nTSTEPS);
xAug(:,1) = [x0; 0]; % cumulative objective value becomes last state

tHalf = t(1:end-1) + h/2;
uHalf = .5*(u(:,1:end-1) + u(:,2:end));

   function value = fullRHS(tVec, xAugVec, uVec)
      xVec = xAugVec(1:end-1,:);
      value = [prob.stateRHS(tVec, xVec, uVec); 
               prob.objective(tVec, xVec, uVec)];
   end

for i = 1:nTSTEPS-1
   k1 = fullRHS(t(i), xAug(:,i), u(:,i));
   k2 = fullRHS(tHalf(i), xAug(:,i) + .5*h*k1, uHalf(:,i));
   k3 = fullRHS(tHalf(i), xAug(:,i) + .5*h*k2, uHalf(:,i));
   k4 = fullRHS(t(i+1), xAug(:,i) + h*k3, uHalf(:,i+1));
   xAug(:,i+1) = xAug(:,i) + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end

x = xAug(1:end-1, :);
J = xAug(end,end);

end

