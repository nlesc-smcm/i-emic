system('tail -n +2 cdata.txt > cdata.tmp'); 
cdata = load('cdata.tmp')

for i = 2:size(cdata,2)
    figure(i)
    plot(cdata(:,1),cdata(:,i));
    title(i)
    xlabel('par');
    grid on;

end
