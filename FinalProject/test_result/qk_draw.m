load('qk_1_1_vector.mat');
load('qk_1_1_vector_new.mat');
load('qk_3_3_vector.mat');
load('qk_3_3_vector_new.mat');
load('qk_5_5_vector.mat');
load('qk_5_5_vector_new.mat');
qk_1_1_vector(1,1)=floor(qk_1_1_vector(1,1)/20);
qk_1_1_vector=qk_1_1_vector./sum(qk_1_1_vector);
qk_1_1_vector_new(1,1)=floor(qk_1_1_vector_new(1,1)/20);
qk_1_1_vector_new=qk_1_1_vector_new./sum(qk_1_1_vector_new);
qk_3_3_vector(1,1)=floor(qk_3_3_vector(1,1)/20);
qk_3_3_vector=qk_3_3_vector./sum(qk_3_3_vector);
qk_3_3_vector_new(1,1)=floor(qk_3_3_vector_new(1,1)/20);
qk_3_3_vector_new=qk_3_3_vector_new./sum(qk_3_3_vector_new);
qk_5_5_vector(1,1)=floor(qk_5_5_vector(1,1)/20);
qk_5_5_vector=qk_5_5_vector./sum(qk_5_5_vector);
qk_5_5_vector_new(1,1)=floor(qk_5_5_vector_new(1,1)/20);
qk_5_5_vector_new=qk_5_5_vector_new./sum(qk_5_5_vector_new);

figure1 = figure('Color',[1 1 1]);
subplot(1,2,1);
hold on 

axis([0 20 0 0.4]);

plot(qk_1_1_vector,'og');
plot(qk_3_3_vector,'or','MarkerFaceColor','r');
plot(qk_5_5_vector,'^b');
legend('(1,1)', '(3,3)', '(5,5)');
for i=1:20
    y=linspace(0,qk_1_1_vector(1,i));
    plot(i,y,'-.g','LineWidth',0.5);
end

for i=1:20
    y=linspace(0,qk_3_3_vector(1,i));
    plot(i,y,'--r');
end

for i=1:20
    y=linspace(0,qk_5_5_vector(1,i));
    plot(i,y,'-b');
end

set(gca, 'XTick', [0 4 8 12 16 20]) ;
set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'});

xlabel('P(w_s|X)');
ylabel('q[k]');
title('Brian Method');

subplot(1,2,2);
hold on 
axis([0 20 0 0.4]);

plot(qk_1_1_vector_new,'og');
plot(qk_3_3_vector_new,'or','MarkerFaceColor','r');
plot(qk_5_5_vector_new,'^b');
legend('(1,1)', '(3,3)', '(5,5)');
for i=1:20
    y=linspace(0,qk_1_1_vector_new(1,i));
    plot(i,y,'-.g','LineWidth',0.5);
end

for i=1:20
    y=linspace(0,qk_3_3_vector_new(1,i));
    plot(i,y,'--r');
end

for i=1:20
    y=linspace(0,qk_5_5_vector_new(1,i));
    plot(i,y,'-b');
end

set(gca, 'XTick', [0 4 8 12 16 20]) ;
set(gca,'XTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0'});

xlabel('P(w_s|X)');
ylabel('q[k]');
title('Proposed');
