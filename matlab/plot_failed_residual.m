[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');
surfm      = landm(2:n+1,2:m+1,l+1);  %Only interior surface points
resa = load('failed_rhs.atmos');
reso = load('failed_rhs.ocean');

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
mx = max(max(abs(To)));
[j,i] = find(abs(To) == mx);
fprintf('maximum in To at i=%d, j=%d\n', i, m+1-j);

figure(1); imagesc(To); title('surface temperature residual'); colorbar
hold on; contour(mx*LM,2); hold off
figure(2); imagesc(TL); title(['temperature residual in layer ', num2str(L)]); colorbar
mx = max(max(abs(TL)));
hold on; contour(mx*LM,2); hold off
figure(3); imagesc(So); title('surface salinity residual');  colorbar
mx = max(max(abs(So)));
hold on; contour(mx*LM,2); hold off
figure(4); imagesc(Ta); title('atmosphere temperature residual');colorbar
mx = max(max(abs(Ta)));
hold on; contour(mx*LM,2); hold off
