function g = test2(F)
global r
g = @gFunc;

   function value = gFunc(x)
      value = F(x,r);
   end


end

