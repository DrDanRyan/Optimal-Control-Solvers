function output = unpack(vector)
   output = cell(1, length(vector));
   for i = 1:length(vector)
      output{i} = vector(i);
   end
end

