figure1 = figure('Color',[1 1 1]);

hold on 
axis([0 20 0 0.5]);

y=linspace(0,0.5);

plot(1,y,'r');
plot(1,0.5,'or');
plot(20,y,'r');
plot(20,0.5,'or');

for i=2:19
    plot(i,0,'or');
end

set(gca, 'XTick', [0 4 8 12 16 20]) ;
set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'});

xlabel('P(w_s|X)');
ylabel('q[k]');
