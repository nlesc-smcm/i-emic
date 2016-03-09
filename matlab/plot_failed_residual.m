[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = readfort44('fort.44');

resa = importdata('failed_rhs.atmos');
reso = importdata('failed_rhs.ocean');
reso = reso(:,3);
norm(reso)
norm(resa)
norm([reso;resa]');

To = zeros(m,n);
Ta = zeros(m,n);

dim = m*n*l*nun;

for j = 1:m
  for i = 1:n
	To(m-j+1,i) = reso(find_row(nun,n,m,l,i,j,l,5));
	Ta(m-j+1,i) = resa(i + (j-1)*m);
  end
end

subplot(2,1,1); imagesc(To); colorbar
subplot(2,1,2); imagesc(Ta); colorbar

