load('mu_det_vector_new.mat');
load('mu_det_vector.mat');

figure1 = figure('Color',[1 1 1]);
hold on 
axis([0 15 0.5 1]);
plot(mu_det_vector, '-ob','LineWidth',2);
plot(mu_det_vector_new,'-^r','LineWidth',2);
legend('Brian Method', 'Proposed');
y=linspace( 0.5, 1);
plot(1,y,'--k');
plot(6,y,'--k');
plot(10,y,'--k');
plot(13,y,'--k');
plot(15,y,'--k');
set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]) ;
set(gca,'XTickLabel',{'(1,1)','(1,2)','(1,3)','(1,4)','(1,5)','(2,2)', ...
    '(2,3)','(2,4)','(2,5)','(3,3)','(3,4)','(3,5)','(4,4)','(4,5)','(5,5)'});
ylabel ('u_{det}');
pos=axis;
xlabel('Segment Size Indices','position',[0.5*pos(2) 0.45]);
rotateticklabel(gca, 90);