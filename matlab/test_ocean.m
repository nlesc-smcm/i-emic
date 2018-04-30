[n m l la nun xmin xmax ymin ymax hdim x y z xu yv zw landm] = ...
    readfort44('fort.44');
dim = n*m*l*nun;

idx = [];
for i = 1:nun;
    idx = [idx, i:nun:dim];
end

idxu = [];
for i = 1:nun;
    idxu(i,:) = [i:nun:dim];
end

jjname = 'ocean_jac';
%jjname = 'mat';
%mmname = 'ocean_B';
icname = 'intcond_coeff';
%ints   = 'integrals2';

%B    = load(mmname);
IC   = load(icname);
%ints = load(ints); 
C    = load(jjname); 
C    = spconvert(C);

IC(22)

Cr = C(idx,idx); % reordering


XX = 6;
Ar = C(idxu(XX,:),idxu(XX,:));
%Ar = Ar(1:end-1,:) % cutoff final integral row
IC = IC(idxu(XX,:)); 
%IC = IC(1:end-1); % cutoff final row

% figure(1); 
% spy(C);
% nnz(C);
% 
% figure(2); 
% spy(Cr);
% 
% figure(3);
% spy(Ar);

intAr = sum(diag(IC(1:end-1))*Ar(1:end-1,:),1);

intAr = reshape(full(intAr), n, m, l);
cma = max(intAr(:));
cmi = min(intAr(:));
for k=1:l
    figure(3+k)
    imagesc(intAr(:,:,k)');
    colorbar;
    title(['level ', num2str(k)]);
    set(gca,'ydir', 'normal')
    % caxis([cmi, cma])
    drawnow
end

% jjname = 'ocean_jac_cont_-';
% mmname = 'ocean_B_cont_-';

% jjname2 = 'ocean_jac_cont_+';
% mmname2 = 'ocean_B_cont_+';

% C = load(jjname); 
% C = spconvert(C);

% C2 = load(jjname2); 
% C2 = spconvert(C2);

% %JnC = load_numjac('ocean_numjac');

% B  = load(mmname);
% B2 = load(mmname2);



% %B  = -ones(1, dim);
% B2 = B;
% B2(B2 == 0) = -1e-1;

% B  = spdiags(B', 0, dim, dim);
% B2 = spdiags(B2', 0, dim, dim);

% idx = [];

% for i = 1:nun;
%     idx = [idx, i:nun:dim];
% end

% Cr   = C(idx,idx); % reordering
% Cr2  = C2(idx,idx);

% figure(1); 
% spy(Cr);

% figure(2);
% Df=(Cr2-Cr)./Cr2 > 1e-2;
% spy(Df);

% figure(3); 
% tic
% opts.maxit = 1000;
% opts.tol = 1e-12;

% neig = 100;
% V = []; D = [];
% [V,D]=eigs(C,B2,neig,0,opts);
% toc

% for i = 1:neig
%     fprintf('%3d: real: %2.5e imag: %2.5e \n', i, real(D(i,i)), ...
%             imag(D(i,i)));
% end

% plot(real(diag(D)),imag(diag(D)),'*');
% axis auto
% grid on;
% hold on;
% plot(xlim, [0,0], 'k-', 'linewidth', 1.1)
% plot([0,0], ylim, 'k-', 'linewidth', 1.1);

% tic
% V2 = []; D2 = [];
% [V2,D2]=eigs(C2,B2,neig,0,opts);
% toc

% for i = 1:neig
%     fprintf('%3d: real: %2.5e imag: %2.5e \n', i, real(D2(i,i)), ...
%             imag(D2(i,i)));
% end

% plot(real(diag(D2)),imag(diag(D2)),'o');
% axis auto
% grid on;
% hold on;
% plot(xlim, [0,0], 'k-', 'linewidth', 1.1)
% plot([0,0], ylim, 'k-', 'linewidth', 1.1);
% hold off;