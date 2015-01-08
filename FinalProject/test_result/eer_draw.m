load('eer_vector_new.mat');
load('eer_vector_orig.mat');
%x=linspace(0,1,15);
figure1 = figure('Color',[1 1 1]);
hold on 
axis([0 15 0 0.5]);
plot(eer_vector_orig, '-ob','LineWidth',2);
plot(eer_vector_new,'-^r','LineWidth',2);
y=linspace(0,0.5);
plot(1,y,'--k');
plot(6,y,'--k');
plot(10,y,'--k');
plot(13,y,'--k');
plot(15,y,'--k');
set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]) ;
set(gca,'XTickLabel',{'(1,1)','(1,2)','(1,3)','(1,4)','(1,5)','(2,2)', ...
    '(2,3)','(2,4)','(2,5)','(3,3)','(3,4)','(3,5)','(4,4)','(4,5)','(5,5)'});
ylabel('EER');
pos=axis;
xlabel('Segment Size Indices','position',[0.5*pos(2) -0.045]);
legend('Brian Method', 'Proposed');

rotateticklabel(gca, 90);
