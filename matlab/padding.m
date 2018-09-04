function [mat] = padding(mat, dim1, dim2)

    [n,m] = size(mat);

    if (n < dim1)
        mat = [mat ; zeros(dim1-n, m)];
    end
    
    if (m < dim2)
        mat = [mat, zeros(dim1, dim2 - m)];
    end        
end