load_atmosphere

% testing matrix vector product with arrays beg, ico and jco

%## C code
%## int first = beg_[row-1];
%## int last  = beg_[row] - 1;
%## double result = 0.0;
%## for (int j = first - 1; j <= last; ++j)
%##   result += ico_[j] * (*state_)[jco_[j]-1];
%## return result;
N = 16*16;
result = zeros(N,1);
for row = 1:N
  first  = beg(row);
  last   = beg(row+1)-1;
  for j = first:last
	result(row) = result(row) + ico(j)*state(jco(j));
  end
end

imagesc(reshape(result,16,16));
