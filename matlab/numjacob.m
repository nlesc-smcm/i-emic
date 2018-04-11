function [J] = numjacob(fun,u,varargin)
N       = length(u);
epsilon = 1e-7;
v       = zeros(N,1);
J       = zeros(N);
f       = feval(fun, u, varargin{:});
for j = 1:N
    v(j)   = 1;
    J(:,j) = (feval(fun, u + epsilon*v, varargin{:}) - f) / (epsilon);
    v(j)   = 0;
end
J = sparse(J); 
