function flipped_array = flip(input_array)
% Inverts array along second dimension. Useful when solving ode backwards in
% time and then interpolating solution vector.

   flipped_array = input_array(:, end:-1:1);
end