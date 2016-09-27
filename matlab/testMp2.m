Mp1 = load('Mp1.ascii'); Mp1 = spconvert(Mp1);
Mp2 = load('Mp2.ascii'); Mp2 = spconvert(Mp2);
Gw  = load('Gw.ascii');  Gw  = spconvert(Gw);
Mp  = load('Mp.ascii');  Mp  = spconvert(Mp);

Gw*Mp'

[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
readfort44('fort.44'); 

irange = 3:6:size(Gw,1);
jrange = 4:6:size(Mp2,2);

Mp  = Mp(:,jrange);

spy(Mp)

Mp2 = Mp(:,end-(n*m)+1:end);
Mp1 = Mp(:,1:end-(n*m));

jd  = logical(max(abs(Mp2),[],1) ~= 0);
figure(1); imagesc(reshape(jd(:),n,m)'); set(gca,'ydir','normal');

size(Mp1)
mMp1 = zeros(n*m,1);
mMp1(jd) = max(abs(Mp1),[],2);

mMp2 = zeros(n*m,1);
mMp2(jd) = max(abs(Mp2),[],2);

figure(2); imagesc(reshape(mMp2,n,m)'); set(gca,'ydir','normal');
title('max(Mp2)')
colorbar
figure(3); imagesc(reshape(mMp1,n,m)'); set(gca,'ydir','normal');
title('max(Mp1)')
colorbar
