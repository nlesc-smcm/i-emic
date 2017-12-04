system('tail -n +2 cdata.txt > cdata.tmp'); 
cdata = load('cdata.tmp');
system('head -n 1 cdata.txt > cdata.tmp');
fid = fopen('cdata.tmp','r');

ncol = size(cdata,2);

titles = cell(ncol,1);

for i = 1:ncol
    titles{i} = fscanf(fid,'%s',1);
end

for i = 2:size(cdata,2)
    figure(i)
    plot(cdata(:,1),cdata(:,i));
    title(titles{i})
    xlabel('par');
    grid on;
end
