function out = array2struct(nms,val)

for jj = 1:length(nms)
   eval(['out.' nms{jj} '=val(jj);']); 
end