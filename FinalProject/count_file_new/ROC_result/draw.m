figure1 = figure('Color',[1 1 1]);
hold on 
plot(diff_vector_3g, same_vector_3g, 'r','LineWidth',2);
plot(diff_vector_2g, same_vector_2g, 'b','LineWidth',2);
plot(diff_vector_4g, same_vector_4g, 'g','LineWidth',2);
plot(diff_vector_5g, same_vector_5g, 'y','LineWidth',2);
plot(diff_vector_6g, same_vector_6g, 'k','LineWidth',2);

x=0:0.1:1;

plot(x,x, '--k');

legend('Proposed', '2-d GMM', '3-d GMM', '4-d GMM', '5-d GMM','Location', 'SouthEast');
xlabel('False positive rate');
ylabel('True positive rate');
