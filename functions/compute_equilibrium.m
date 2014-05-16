function [xStar, lamStar, uStar] = ...
   compute_equilibrium(prob, xGuess, lamGuess, uGuess, r)

   nStates = length(xGuess); % does not include xJ state
   nControls = length(uGuess);
   xIdx = 1:nStates;
   lamIdx = nStates+1:2*nStates;
   uIdx = 2*nStates+1:2*nStates+nControls;
   
   function value = equilibrium_system(y)
      x = y(xIdx); lam = y(lamIdx); u = y(uIdx);
      value = zeros(size(y));
      
      F = prob.F(0, x, u);
      value(xIdx) = F(xIdx);
      
      dFdx_times_vec = prob.dFdx_times_vec(0, x, u, [lam; 1]);
      value(lamIdx) = r*lam - dFdx_times_vec(xIdx);
      
      value(uIdx) = prob.dFdu_times_vec(0, x, u, [lam; 1]);
   end

   yStar = fsolve(@equilibrium_system, [xGuess; lamGuess; uGuess]);
   
   xStar = yStar(xIdx);
   lamStar = yStar(lamIdx);
   uStar = yStar(uIdx);

end

