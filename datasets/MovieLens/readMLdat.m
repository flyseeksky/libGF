function [m_train,m_test,c_userInfo,c_movieInfo] = readMLdat
pth=which('readMLdat');
path_left=pth(1:end-(length('readMLdat')+2));
file=strcat(path_left,'ml-100k\u1');
D = readtable(file,'Delimiter','\t','ReadVariableNames',false);
TrainData=table2array(D);
file=strcat(path_left,'ml-100k\u1_test');
D = readtable(file,'Delimiter','\t','ReadVariableNames',false);
TestData=table2array(D);
file=strcat(path_left,'ml-100k\user');
D = readtable(file,'Delimiter','|','ReadVariableNames',false);
c_userInfo=table2cell(D);
file=strcat(path_left,'ml-100k\u.item.txt');
D = readtable(file,'Delimiter','|','ReadVariableNames',false);
c_movieInfo=table2cell(D);
m_train=zeros(size(c_userInfo,1),size(c_movieInfo,1));
for(i=1:size(TrainData,1))
m_train(TrainData(i,1),TrainData(i,2))=TrainData(i,3);
end

m_test=zeros(size(c_userInfo,1),size(c_movieInfo,1));
for(i=1:size(TestData,1))
m_test(TestData(i,1),TestData(i,2))=TestData(i,3);
end
end

