function fInterp = vectorInterpolant(x, v, interpType)
   nCOMPONENTS = size(v, 1);
   if nCOMPONENTS == 1
      fInterp = griddedInterpolant(x, v, interpType);
   else
      interpArray = cell(nCOMPONENTS, 1);
      for idx = 1:nCOMPONENTS
         interpArray{idx} = griddedInterpolant(x, v(idx, :), interpType);
      end
      fInterp = @(t) cell2mat(cellfun(@(f) f(t), interpArray, 'UniformOutput', false));
   end
end