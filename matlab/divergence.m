function [DIV,UHX,VHY] = divergence(u,v,varz,vary,varx,periodic)

  if nargin < 6
	periodic = false;
  end
  
  %% depth-integration
  for k = 1:length(varz)-1
	u(:,:,k)  = u(:,:,k) *(varz(k+1)-varz(k));
	v(:,:,k)  = v(:,:,k) *(varz(k+1)-varz(k));
  end
  
  UINTH  = sum(u,3);
  VINTH  = sum(v,3);
  UHX    = UINTH;
  VHY    = VINTH;
  
  for j = 1:size(varx,2)
	for i = 1:size(varx,1)-1
	  % assuming periodic bdys
	  if i+1 <= length(varx)-1
		UHX(i,j) = (UINTH(i+1,j) - UINTH(i,j)) / (varx(i+1,j)-varx(i,j));
	  elseif periodic
		UHX(i,j) = (UINTH(1,j) - UINTH(i,j)) / (varx(i+1,j)-varx(i,j));
	  else
		UHX(i,j) = ( -UINTH(i,j)) / (varx(i+1,j)-varx(i,j));
	  end
	end
  end

  for j = 1:length(vary)-1
	if j+1 <= length(vary)-1
	  VHY(:,j) = (VINTH(:,j+1) - VINTH(:,j)) / (vary(j+1)-vary(j));
	else
	  VHY(:,j) = (- VINTH(:,j)) / (vary(j+1)-vary(j));
	end	
  end

  DIV = UHX+VHY;  
end
