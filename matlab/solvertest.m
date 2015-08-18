% make sure that we are in a run directory and that all the outputs are ON
% creating the ocean-atmosphere coupled system

% load / build system
fprintf('loading...\n');

A = load_ocean_matrix;
D = load_atmosphere_matrix;
N = size(A,1);
M = size(D,1);

id = importdata('coupled_rowsB.txt');
vB = importdata('coupled_B.txt');
vC = importdata('coupled_C.txt');

B = sparse(id, 1:M, vB, N, M);
C = sparse(1:M, id, vC, M, N);

% Make A non-singular
A = A + speye(N);

J = [A,B;C,D];


fprintf('\nloading... done\n');

% solver pars:
restart = []; tol = 1e-8; maxit = N+M;

x_exact   =  ones(N+M,1);
x1        =  x_exact(1:N);
x2        =  x_exact(N+1:N+M);
b         =  J * x_exact;
b1        =  b(1:N);
b2        =  b(N+1:N+M);

% Need to think of something iterative and magical....
X1 = rand(N,1); %seed
norm(X1 - x1)
X2 = D \ (b2 - C*X1);
norm(X2 - x2)
for i = 1:10
X2 = D \ (b2 - C*inv(A)*B*X2);
norm(X2 - x2)
end
X1 = A \ (b1 - B*X2); % I want just 1 solve with A and an accurate solution...
norm(X1 - x1)




%## % ELIMINATION SOLVE
%## z1  = gmres(A,b1,restart,tol,maxit);
%## c2  = b2 - C*z1;
%## z2  = gmres(D,c2,restart,tol,maxit);
%## z3  = gmres(A,B*z2,restart,tol,maxit);
%## z4  = gmres(D,C*z3,restart,tol,maxit);
%
%## X2 = z2 + z4;
%
%## X1  = gmres(A, (b1-B*X2), restart, tol, maxit);

