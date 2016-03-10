[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
surfm      = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
%resa = importdata('failed_rhs.atmos');
%reso = importdata('failed_rhs.ocean');
%reso = reso(:,3);

norm(reso)
norm(resa)
norm([reso;resa]')

L = l-5;

To = zeros(m,n);
TL = zeros(m,n);
So = zeros(m,n);
Ta = zeros(m,n);

dim = m*n*l*nun;

for j = 1:m
  for i = 1:n
	To(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,l,5));	
	TL(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,L,5));
	LM(m-j+1,i) = surfm(i,j);
	So(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,l,6));
	Ta(m-j+1,i) = resa(i + (j-1)*m);
  end
end

[j,i] = find(To==max(max(To)));
fprintf('maximum in To at i=%d, j=%d\n', i, m+1-j);

subplot(4,1,1); imagesc(To); title('surface temperature residual'); colorbar
hold on; contour(1e-4*LM,5); hold off
subplot(4,1,2); imagesc(TL); title(['temperature residual in layer ', num2str(L)]); colorbar
subplot(4,1,3); imagesc(So); title('surface salinity residual');  colorbar
subplot(4,1,4); imagesc(Ta); title('atmosphere temperature residual');colorbar
