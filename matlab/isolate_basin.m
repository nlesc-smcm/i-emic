function [new_mask] = isolate_basin(original_mask, level)
		 
  maskp = original_mask;
  [m,n,l] = size(maskp);

  if (l > 1)
	fprintf('input mask is not in single layer format!');
	return;
  end
  
  masktmp = maskp;
  new_mask = maskp;

  if nargin < 2
	level = 1; % division level default=1
  end
  
  masktmp(masktmp < level) = 0;
  masktmp(masktmp >= level) = 1;

  % double mask to incorporate periodicity
  maskperiod = [masktmp(:,round(n/2):end),masktmp,masktmp(:,1:round(n/2))];
  id  = numel(round(n/2):n)+1:numel(round(n/2):n)+n;

  % obtain connected components
  CC = bwconncomp(maskperiod,4);

  % this is tricky: we pick the first idxlist because up until now this
  % has always been the right basin.
  id1 = CC.PixelIdxList{1}; 
  
  maskperiod(id1) = 999;  
  new_mask(maskperiod(:,id)~=999) = 0;  
end
