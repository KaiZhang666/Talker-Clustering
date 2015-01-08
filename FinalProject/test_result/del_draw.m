load('del_vector_new.mat');
load('del_vector.mat');

del=zeros(1,15);
del_new=zeros(1,15);
for i=1:15
    del_sum=sum(del_vector(1,(i-1)*40+21:(i-1)*40+40));
    del_sum_new=sum(del_vector_new(1,(i-1)*40+21:(i-1)*40+40));
    for j=1:20
        del_vector((i-1)*40+j)=del_vector((i-1)*40+j)/del_vector((i-1)*40+j+20);
        del_vector_new((i-1)*40+j)=del_vector_new((i-1)*40+j)/del_vector_new((i-1)*40+j+20);
        del_vector((i-1)*40+j+20)=del_vector((i-1)*40+j+20)/del_sum;
        del_vector_new((i-1)*40+j+20)=del_vector_new((i-1)*40+j+20)/del_sum_new;      
    end
end
for i=1:15
    for j=1:20
        del(1,i)=del(1,i)+del_vector((i-1)*40+j+20)*(del_vector((i-1)*40+j)-0.025-(j-1)*0.05)^2;
         del_new(1,i)=del_new(1,i)+del_vector_new((i-1)*40+j+20)*(del_vector_new((i-1)*40+j)-0.025-(j-1)*0.05)^2;
    end
end
       
figure1 = figure('Color',[1 1 1]);
hold on 
axis([0 15 0 0.1]);
plot(del, '-ob','LineWidth',2);
plot(del_new,'-^r','LineWidth',2);
legend('Brian Method', 'Proposed');
y=linspace(0,0.1);
plot(1,y,'--k');
plot(6,y,'--k');
plot(10,y,'--k');
plot(13,y,'--k');
plot(15,y,'--k');
set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]) ;
set(gca,'XTickLabel',{'(1,1)','(1,2)','(1,3)','(1,4)','(1,5)','(2,2)', ...
    '(2,3)','(2,4)','(2,5)','(3,3)','(3,4)','(3,5)','(4,4)','(4,5)','(5,5)'});
ylabel ('\delta_{err}');
pos=axis;
xlabel('Segment Size Indices','position',[0.5*pos(2) -0.01]);
rotateticklabel(gca, 90);    
    
    
