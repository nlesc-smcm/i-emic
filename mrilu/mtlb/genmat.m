% Generate the sparse square matrices  A  and  B, a lower triangular
% matrix  L  and a vector D containing the nonzero diagonal of a diagonal
% matrix and a vector  x.
% Write the matrices in CSR format to the files  matA  and  matB.
% Write the vector  x  to the file  vecx
% Uses the following global variables
%  n        the order of the matrices and vectors
%  density  the percentage of the nonzero elements in the matrices
%
A = sprandn(n, n, density);
[begA,jcoA,coA] = sp2csr(A);
wrtbcsr(begA,jcoA,coA, 'matA');
% 
B = sprandn(n, n, density);
[begB,jcoB,coB] = sp2csr(B);
wrtbcsr(begB,jcoB,coB, 'matB');
%
x = rand(n,1);
wrtbvec(x, 'vecx');
%
C = sprandsym(n, density, 1+rand(1,n));
U = chol(C);
D = full(diag(U));
R = spdiags(1./D,0,n,n) * U;
D = D.^2;
wrtbvec(D, 'diagD');
%
L = R';
%
LDLt = L*spdiags(D,0,n,n)*L';
% 
Rest = norm(LDLt - C, inf);
%
L = L - spdiags(full(diag(L)),0,n,n);
[begL,jcoL,coL] = sp2csr(L);
wrtbcsr(begL,jcoL,coL, 'matL');
