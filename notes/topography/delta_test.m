N     = 10;
M     = 100;
delta = linspace(0,N-1,M);

b  = zeros(N,1);
a  = d;
kr = d;

a(1) = 0;
b(1) = 1;

for i = 1:M-1
  k = floor(delta(i));
  kr(k+1) = k;
  a(k+2) = mod(k + 1,2)*2 + a(k+1);
  b(k+2) = mod(k, 2)*2 + b(k+1);
end

[kr,a,b]
