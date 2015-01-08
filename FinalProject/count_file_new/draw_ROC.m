sameMatrix_6g=[];
diffMatrix_6g=[];

same_list=dir('D:\Thesis\FinalProject\count_file_new\same_6g\*.txt');
[length,width]=size(same_list);
for i=1:length
    same_listname=same_list(i).name;
    sameMatrix_6g=[sameMatrix_6g; load(char(['D:\Thesis\FinalProject\count_file_new\same_6g\' same_listname]))];
end

diff_list=dir('D:\Thesis\FinalProject\count_file_new\diff_6g\*.txt');
[length,width]=size(diff_list);
for i=1:length
    diff_listname=diff_list(i).name;
    diffMatrix_6g=[diffMatrix_6g; load(char(['D:\Thesis\FinalProject\count_file_new\diff_6g\' diff_listname]))];
end

sameMatrix_6g=sum(sameMatrix_6g);
diffMatrix_6g=sum(diffMatrix_6g);

same_vector_6g=zeros(1,100);
diff_vector_6g=zeros(1,100);
for i=1:100
    same_vector_6g(1,i)=sameMatrix_6g(1,(300-3*i+1))/sameMatrix_6g(1,(300-3*i+3));
    diff_vector_6g(1,i)=diffMatrix_6g(1,(300-3*i+2))/diffMatrix_6g(1,(300-3*i+3));
end
for i=1:99
    if same_vector_6g(1,i+1)<same_vector_6g(1,i)
        same_vector_6g(1,i+1)=same_vector_6g(1,i);
    end
    if diff_vector_6g(1,i+1)<diff_vector_6g(1,i)
        diff_vector_6g(1,i+1)=diff_vector_6g(1,i);
    end   
end
same_vector_6g(100)=1; diff_vector_6g(100)=1;
line=0:0.01:1;
hold on
plot(line,line,'--k');
plot(diff_vector_6g, same_vector_6g,'b');
    