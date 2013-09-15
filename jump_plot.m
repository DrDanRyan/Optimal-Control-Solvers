function jump_plot(jumpPtArray, yValues, varargin)
% Plots the piecewise constant function with n segments defined by the row
% vector jumpPtArray (lenght n+1) and the column vector yValues (length n).
% The third argument is an optional color string.

if size(varargin) > 0
   color = varargin{1};
else
   color = 'b';
end

plot([jumpPtArray(1:end-1); jumpPtArray(2:end)], [yValues'; yValues'], color);