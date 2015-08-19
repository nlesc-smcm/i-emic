function [X1,X2] = solvertest(relaxation);
  % In this function we test a backward block GS solver with overrelaxation
  % for the coupled ocean+atmosphere problem
		 
  % make sure that we are in a run directory and that all the outputs are ON
		 
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

  % Make A diagonally dominant non-singular
  A = A + 1*speye(N);;

  J = [A,B;C,D];

  fprintf('\nloading... done\n');

  % solver pars:
  w = relaxation;

  x_exact   =  ones(N+M,1);
  xx1       =  x_exact(1:N);
  xx2       =  x_exact(N+1:N+M);
  b         =  J * x_exact;
  b1        =  b(1:N);
  b2        =  b(N+1:N+M);

  X1 = rand(N,1); %fprintf('norm(X1 - xx1): %e\n', norm(X1 - xx1));
  X2 = rand(M,1); %fprintf('norm(X2 - xx2): %e\n', norm(X2 - xx2));

  % block backward SOR
  X2 = D \ (-w*C*X1 + (1-w)*D*X2 + w*b2);
  fprintf('1: %e\n', norm(X2 - xx2));

  X1 = A \ (-w*B*X2 + (1-w)*A*X1 + w*b1);
  fprintf('2: %e\n', norm(X1 - xx1));

  X2 = D \ (-w*C*X1 + (1-w)*D*X2 + w*b2);
  fprintf('3: %e\n', norm(X2 - xx2));
  
end
