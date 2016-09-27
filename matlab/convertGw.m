% Gw and Mp matrices
Gw = load('Gw.ascii');
Gw = spconvert(Gw);
Mp = load('Mp.ascii');
Mp = spconvert(Mp);

% convert to correct data
[MGw,NGw] = size(Gw);
[MMp,NMp] = size(Mp);

m = 16; n = 16; l = 8;

id  = logical(max(abs(Gw),[],2) ~= 0);
jd  = logical(max(abs(Mp),[],1) ~= 0);
Gw  = Gw(id,jd);
Mp  = Mp(:,jd);
Gw1 = Gw(:,1:nnz(id));
Gw2 = Gw(:,nnz(id)+1:end);
Mp1 = Mp(:,1:nnz(id));
Mp2 = Mp(:,nnz(id)+1:end);

fprintf('Gw1*Mp1 + Gw2*Mp2 = %f\n',...
		norm(Gw1*Mp1' + Gw2*Mp2',Inf));

% get b
b    = load('b.ascii');
bhat = load('bhat.ascii');
w    = load('w.ascii');
u    = load('u.ascii');
v    = load('v.ascii');
z    = load('z.ascii');
wmz  = load('wmz.ascii');
x    = load('x.ascii');

W = Gw1 \ b';
fprintf('w-W     = %f\n', norm(w'-W));

U = Mp1*W;
fprintf('u-U     = %f\n', norm(u'-U));

Z = Mp1'*U;
fprintf('z-Z     = %f\n', norm(z'-Z));

WMZ = W-Z;
fprintf('wmz-WMZ = %f\n', norm(wmz'-WMZ));

V = -Mp2'*U;
fprintf('v-V     = %f\n', norm(v'-V));

X = [WMZ;V];
fprintf('x-X     = %f\n', norm(x'-X));

format long
norm(Gw1*WMZ)
norm(Gw2*V)

norm(b' - Gw1*WMZ - Gw2*V)
norm(Mp1*WMZ+Mp2*V)

%norm([b';zeros(nnz(jd)-nnz(id),1)]-[Gw;Mp]*x')

norm(WMZ)
norm(V)
norm([WMZ;V])
