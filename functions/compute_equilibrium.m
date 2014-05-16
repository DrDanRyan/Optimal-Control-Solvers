function [xStar, lamStar, uStar] = ...
   compute_equilibrium(prob, xGuess, lamGuess, uGuess)

   nStates = length(xGuess);
   nControls = length(uGuess);
   
   function value = equilibrium_system(y)
      x = y(1:nStates); lam = y(nStates+1:2*nStates); u = y(2*nStates+1:end);
      value = zeros(2*nStates+nControls, 1);
      F = prob.F(0, x, u);
      value(1:nStates) = F(1:end-1,:);
      dFdx_times_vec = prob.dFdx_times_vec(0, x, u, lam);
      value(nStates+1:2*nStates) = dFdx_times_vec(1:end-1,:);
      value(2*nStates+1:end) = prob.dFdu_times_vec(0, x, u, lam);
   end

   yStar = fsolve(@equilibrium_system, [xGuess; lamGuess; uGuess]);
   
   xStar = yStar(1:nStates);
   lamStar = yStar(nStates+1:2*nStates);
   uStar = yStar(2*nStates+1:end);

end

