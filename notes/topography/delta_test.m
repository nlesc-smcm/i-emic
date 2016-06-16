N     = 10;
M     = 100;
delta = linspace(0,N-1,M);

d = zeros(N,1);
s = d;
d(1) = 1;
s(1) = 0;
for i = 1:M-1
  k = floor(delta(i));
  d(k+2) = mod(floor(delta(i)), 2)*2 + d(k+1);
  s(k+2) = mod(floor(delta(i) + 1),2)*2 + s(k+1);
end

[d,s]
