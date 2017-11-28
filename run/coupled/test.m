cdata = dlmread('cdata.txt', ' ', 1, 0)

zero_idx = sum(cdata,1) == 0;
cdata = cdata(:,~zero_idx);

for i = 2:size(cdata,2)
    figure(i)
    plot(cdata(:,1),cdata(:,i))
    title(i)
    xlabel('par')
    grid on

end
