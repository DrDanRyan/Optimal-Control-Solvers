function value = heval(func, tspan, components)
% Evaluates the vector valued funciton handle at the points in tspan but
% only returns the components specified in components.
   full_vec = func(tspan);
   value = full_vec(components, :);
end