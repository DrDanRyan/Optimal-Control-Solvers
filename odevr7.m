function [tout,yout] = odevr7(odefun,tspan,y0,re,ae)
% Solve initial value problems y' = f(t,y).  The user interface is much
% like the solvers of the Matlab ODE Suite. **ODEFUN must be coded so that 
% f(T,Y) returns YP with YP(:,k) = f(T(k),Y(:,k)) for each entry in T and 
% corresponding column in Y.** The relative error tolerance RE is a scalar
% that must be positive. It has default  value 1e-3.  If necessary, the
% input RE will be changed so that 100*EPS <= RE <= 0.01. The absolute 
% error tolerance AE can be a scalar or vector. All components of AE must 
% be positive. The default is a scalar with value 1e-6.

t = tspan(1);
tend = tspan(end);
dir = sign(tend - t);
tout_given = length(tspan) > 2;
y = y0(:);
neq = length(y);
yp = odefun(t,y);

if nargin < 4 || isempty(re)
    re = 1e-3;
elseif re <= 0
    error('MATLAB:RENotPos', 'RE must be positive.');
elseif re > 0.01
    re = 0.01;
    warning('MATLAB:REdecrease','RE has been decreased to 0.01.')
end
if re < 100 * eps 
  re = 100 * eps;
  warning('MATLAB:REIncrease','RE has been increased to %g.',re)
end
if nargin < 5 || isempty(ae)
    ae = 1e-6;
elseif any(ae <= 0)
  error('MATLAB:AENotPos', 'AE must be positive.');
end
% Redefine AE for use in forming WT:
ae = ae(:)/re;

% This method of global order 7 computes a block of 7 values in each step. 
pow = 1/8; 
BlockSize = 7;         

% Compute an initial step size h using y'(t).
hmax = 0.1*abs(tend - t);
hmin = 16*max(eps(t),eps(tend));
rh = norm((yp ./ max(abs(y),(ae/re))),inf) / (0.8*re^pow);
absh = hmax;
if absh*rh > 1
    absh = 1/rh;
end
h = dir*max(absh,hmin);

if tout_given
    Lout = length(tspan);
    tout = tspan(:);
    yout = zeros(Lout,neq);
else
    Lout = 100;
    tout = zeros(Lout,1);
    yout = zeros(Lout,neq);
end

% Initialize output.
nout = 1;
tout(nout) = t;
yout(nout,:) = y';

% Define coefficients of method. 
M = [ 139849/846720, 1466/6615, 1359/6272, 1448/6615, ...
      36725/169344, 54/245, 3577/17280 ; ...
     -4511/31360, -71/2940, 1377/31360, 8/245, ...
      775/18816, 27/980, 49/640 ; ...
      123133/846720, 68/735, 5927/31360, 1784/6615, ...
      4625/18816, 68/245, 2989/17280 ; ...
     -88547/846720, -1927/26460, -3033/31360, -106/6615, ...
      13625/169344, 27/980, 2989/17280 ; ...
      1537/31360, 26/735, 1377/31360, 8/245, ...
      1895/18816, 54/245, 49/640 ; ...
     -11351/846720, -29/2940, -373/31360, -64/6615, ...
     -275/18816, 41/980, 3577/17280 ; ...
      275/169344, 8/6615, 9/6272, 8/6615, 275/169344, ...
      0, 751/17280 ];  
CoefT = [1/7, 2/7, 3/7, 4/7, 5/7, 6/7, 1];
Coefyp = [751/17280, 41/980, 265/6272, ...
          278/6615, 265/6272, 41/980, 751/17280];
 
done = false;
while ~done
    
    % Adjust step size near the end of the interval.
    if 1.1*abs(h) >= abs(tend - t)     % stretch
        h = tend - t;
        done = true;
    elseif 2*abs(h) >= abs(tend - t)   % look-ahead
        h = (tend - t)/2;
    end

    failed = false;
    while true  % Reduce h until the step to t + h succeeds.
        absh = abs(h);
        T = t + h * CoefT; 
        temp = y(:,ones(1,BlockSize)) + yp * ( h * Coefyp);   
        hM = h * M;
        Y = temp + yp(:,ones(1,BlockSize)) * hM; 
        for k = 1:6
            YP = odefun(T,Y);
            Y = temp + YP * hM;
        end
             
        % Control local error of approximation at each point of T and the 
        % scaled residual at those points.  YP holds the derivative of the 
        % polynomial solution at points of T. Compute the residual using 
        % the values Y at points of T in the ODEs. Note that if the step 
        % is accepted, the slope at t+h, DERIV(:,end), can be used for the 
        % next step and the computation of the error estimate is "free".  
        DERIV = odefun(T,Y);
        wt = ae + max(abs(Y),[],2);
        errnrm = absh * max( max( abs(DERIV - YP) ./ wt(:,ones(1,BlockSize) ) ) );
        
        if errnrm <= re  % successful step
            break
        else             % reduce h, try again
            if absh < 16*max(eps(t),eps(tend))
                error('Step size needed is too small for the precision.')
            end
            h = h*max(0.5,0.8*(re/errnrm)^pow);
            failed = true;
            done = false;
        end  
        
    end  % Loop until step succeeds.
    
    % If TOUT is given, see whether any output points lie in (t,t+h].
    % If so, evaluate the continuous extension to get solution there.
    if tout_given
        indices = find( (dir*(tout - t) > 0) & (dir*(t+h - tout) >= 0) );
        if ~isempty(indices)
            yint = continuous_extension(tout(indices));
            numout = length(indices);                        
            yout(nout+1:nout+numout,:) = yint';
            nout = nout + numout;
        end
    else
        if nout + BlockSize > Lout       % Allocate more storage as needed.
            Lout = Lout + 100;
            tout(Lout) = 0;
            yout(Lout,1:neq) = 0;
        end
        tout(nout+1:nout+BlockSize) = T';
        yout(nout+1:nout+BlockSize,:) = Y';
        nout = nout + BlockSize;
    end
    
    if ~done                 
        t = T(end);
        y = Y(:,end);
        yp = DERIV(:,end);
        if ~failed           % No increase after a failed attempt.
            absh = absh/max(0.2,1.25*(errnrm/re)^pow);
            h = sign(h)*min(absh,hmax);
        end
    end
      
end    

% Trim output arrays.
tout = tout(1:nout);
yout = yout(1:nout,:);


%===Nested function========================================================
function yce = continuous_extension(tint)
    s = ( tint(:)' - t )/h; 
    s2 = s .* s; s3 = s .* s2; s4 = s .* s3; s5 = s .* s4;
    s6 = s .* s5; s7 = s .* s6; s8 = s .* s7;
    yce = y(:,ones(size(tint))) ...
        - (117649*h/720)*yp*((1/8)*s8-(4/7)*s7+(23/21)*s6-(8/7)*s5+...
        + (967/1372)*s4-(268/1029)*s3+(6534/117649)*s2-(720/117649)*s) ...
        + (823543*h/720)*YP(:,1)*((1/8)*s8-(27/49)*s7+(295/294)*s6-...
          (333/343)*s5+(1276/2401)*s4-(2676/16807)*s3+(360/16807)*s2) ...
        - (823543*h/240)*YP(:,2)*((1/8)*s8-(26/49)*s7+(45/49)*s6-...
          (284/343)*s5+(3929/9604)*s4-(1758/16807)*s3+(180/16807)*s2) ...
        + (823543*h/144)*YP(:,3)*((1/8)*s8-(25/49)*s7+(247/294)*s6-...
          (1219/1715)*s5+(778/2401)*s4-(3796/50421)*s3+(120/16807)*s2) ...
        - (823543*h/144)*YP(:,4)*((1/8)*s8-(24/49)*s7+(113/147)*s6-...
          (1056/1715)*s5+(2545/9604)*s4-(984/16807)*s3+(90/16807)*s2)...
        + (823543*h/240)*YP(:,5)*((1/8)*s8-(23/49)*s7+(69/98)*s6-...
          (185/343)*s5+(536/2401)*s4-(804/16807)*s3+(72/16807)*s2) ...
        - (823543*h/720)*YP(:,6)*((1/8)*s8-(22/49)*s7+(95/147)*s6-...
          (164/343)*s5+(1849/9604)*s4-(2038/50421)*s3+(60/16807)*s2)...
        + (117649*h/720)*YP(:,7)*((1/8)*s8-(3/7)*s7+(25/42)*s6-...
          (3/7)*s5+(58/343)*s4-(12/343)*s3+(360/117649)*s2);
end % continuous_extension
%==========================================================================

end % odevr7
