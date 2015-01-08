eerMatrix=[];
delMatrix=[];
rejMatrix=[];
detMatrix=[];

pdir='D:\Thesis\FinalProject\count_file_new\';

eerlist=dir('D:\Thesis\FinalProject\count_file_new\eer_count*.txt');
dellist=dir('D:\Thesis\FinalProject\count_file_new\del_count*.txt');
rejlist=dir('D:\Thesis\FinalProject\count_file_new\rej_count*.txt');
detlist=dir('D:\Thesis\FinalProject\count_file_new\det_count*.txt');

[length,width]=size(eerlist);
for i=1:length
    eerlistname=eerlist(i).name;
    eerMatrix=[eerMatrix; load(char([pdir eerlistname]))];
end



[length,width]=size(dellist);
for i=1:length
    dellistname=dellist(i).name;
    delMatrix=[delMatrix; load(char([pdir dellistname]))];
end



[length,width]=size(rejlist);
for i=1:length
    rejlistname=rejlist(i).name;
    rejMatrix=[rejMatrix; load(char([pdir rejlistname]))];
end


[length,width]=size(detlist);
for i=1:length
    detlistname=detlist(i).name;
    detMatrix=[detMatrix; load(char([pdir detlistname]))];
end


eer_vector=sum(eerMatrix);
fpr = zeros(1,100);
fnr = zeros(1,100);
eer_offset=4200;
eer_sub_vector=eer_vector(1,eer_offset+1:eer_offset+300);
for i=1 : 3 : 300
    fpr(1,ceil(i/3))=eer_vector(1,i)/(eer_vector(1,i+2));
    if fpr(1,ceil(i/3)) >1
        fpr(1,ceil(i/3))=0;
    end
    fnr(1,ceil(i/3))=eer_vector(1,i+1)/(eer_vector(1,i+2));
    if fnr(1,ceil(i/3)) >1
        fnr(1,ceil(i/3))=0;
    end
end

del_vector=sum(delMatrix);
qk_1_1_vector=del_vector(1,21:40);
Psk_1_1_vector=del_vector(1,1:20);
qk_3_3_vector=del_vector(1,381:400);
Psk_3_3_vector=del_vector(1,361:380);
qk_5_5_vector=del_vector(1,581:600);
Psk_5_5_vector=del_vector(1,561:580);

det_vector=sum(detMatrix);
mu_det_vector=zeros(1,15);
for i=1:15
    mu_det_vector(1,i)= det_vector(1,2*i-1)/det_vector(1,2*i);
end

rej_vector=sum(rejMatrix);
mu_rej_vector=zeros(1,15);
for i=1:15
    mu_rej_vector(1,i)= rej_vector(1,2*i-1)/rej_vector(1,2*i);
end

