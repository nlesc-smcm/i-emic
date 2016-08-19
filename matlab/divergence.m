function [DIV,UHX,VHY] = divergence(u,v,varz,vary,varx)

%% depth-integration
for k = 1:length(varz)-1
    u(:,:,k)  = u(:,:,k)*(varz(k+1)-varz(k));
	v(:,:,k)  = v(:,:,k)*(varz(k+1)-varz(k));
end
UINTH  = sum(u,3);
VINTH  = sum(v,3);
UHX    = UINTH;
VHY    = VINTH;

for i = 1:length(varx)-1
  % assuming periodic bdys
  if i+1 <= length(varx)-1
	UHX(i,:) = (UHX(i+1,:) - UHX(i,:)) / (varx(i+1)-varx(i));
  else
	UHX(i,:) = (UHX(1,:) - UHX(i,:)) / (varx(i+1)-varx(i));
  end
end

for j = 1:length(vary)-1
  if j+1 <= length(vary)-1
	VHY(:,j) = (VHY(:,j+1) - VHY(:,j)) / (vary(j+1)-vary(j));
  else
	VHY(:,j) = (- VHY(:,j)) / (vary(j+1)-vary(j));
  end
	 
end

DIV = UHX+VHY;
