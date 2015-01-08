 load('Psk_1_1_vector.mat');
load('Psk_1_1_vector_new.mat');
load('Psk_3_3_vector.mat');
load('Psk_3_3_vector_new.mat');
load('Psk_5_5_vector.mat');
load('Psk_5_5_vector_new.mat');
 load('qk_1_1_vector.mat');
load('qk_1_1_vector_new.mat');
load('qk_3_3_vector.mat');
load('qk_3_3_vector_new.mat');
load('qk_5_5_vector.mat');
load('qk_5_5_vector_new.mat');

for i=1:20
    Psk_1_1_vector(1,i)=Psk_1_1_vector(1,i)/qk_1_1_vector(1,i);
    Psk_3_3_vector(1,i)=Psk_3_3_vector(1,i)/qk_3_3_vector(1,i);
    Psk_5_5_vector(1,i)=Psk_5_5_vector(1,i)/qk_5_5_vector(1,i);
    Psk_1_1_vector_new(1,i)=Psk_1_1_vector_new(1,i)/qk_1_1_vector_new(1,i);
    Psk_3_3_vector_new(1,i)=Psk_3_3_vector_new(1,i)/qk_3_3_vector_new(1,i);
    Psk_5_5_vector_new(1,i)=Psk_5_5_vector_new(1,i)/qk_5_5_vector_new(1,i);
end


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

Pbin=0.025:0.05:0.975;

figure1 = figure('Color',[1 1 1]);
subplot(1,2,1);
hold on
plot(0.9,0.1,'og');
plot(0.9,0.1,'or','MarkerFaceColor','r');
plot(0.9,0.1,'^b');
legend('(1,1)', '(3,3)', '(5,5)','Location', 'SouthEast' );
for i=1:20
  plot(Pbin(1,i),Psk_1_1_vector(1,i),'og','MarkerSize',100*qk_1_1_vector(1,i));
end
for i=1:20
  plot(Pbin(1,i),Psk_3_3_vector(1,i),'or','MarkerSize',100*qk_3_3_vector(1,i),'MarkerFaceColor','r');
end
for i=1:20
  plot(Pbin(1,i),Psk_5_5_vector(1,i),'^b','MarkerSize',100*qk_5_5_vector(1,i));
end


plot(Pbin,Pbin,'--k');

xlabel('P(w_s|X)');
ylabel('P_s[k]');
title('Brian Method');

subplot(1,2,2);
hold on
plot(0.9,0.1,'og');
plot(0.9,0.1,'or','MarkerFaceColor','r');
plot(0.9,0.1,'^b');
legend('(1,1)', '(3,3)', '(5,5)','Location', 'SouthEast' );
for i=1:20
  plot(Pbin(1,i),Psk_1_1_vector_new(1,i),'og','MarkerSize',100*qk_1_1_vector_new(1,i));
end
for i=1:20
  plot(Pbin(1,i),Psk_3_3_vector_new(1,i),'or','MarkerSize',100*qk_3_3_vector_new(1,i),'MarkerFaceColor','r');
end
for i=1:20
  plot(Pbin(1,i),Psk_5_5_vector_new(1,i),'^b','MarkerSize',100*qk_5_5_vector_new(1,i));
end


plot(Pbin,Pbin,'--k');

xlabel('P(w_s|X)');
ylabel('P_s[k]');
title('Proposed');
