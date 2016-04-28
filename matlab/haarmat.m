function [H] = haarmat(p)
if p == 1
    H = 1;
    return
end
H   = 1/sqrt(2)*[1 1; 1 -1];
dim = 2;
while dim ~= p    
    H = 1/sqrt(2)*[kron(H,[1 1]); kron(eye(dim),[1 -1])];
    dim = size(H,1);
end
H = sparse(H);
end
