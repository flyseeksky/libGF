function [m_clustUser,m_clustUserInfoAge,m_clustUserInfoSex,m_clustUserInfoOccup,m_train,m_test,c_userInfo,c_movieInfo]=prepareMLdat
%% prepare the user info matrix
[m_train,m_test,c_userInfo,c_movieInfo] = readMLdat;
%retain only users that have rated at least 20 movies
ind=1;
for k=1:size(m_train,1)
    if(nnz(m_test(k,:))>20)
        v_usersToKeep(ind)=k;
        ind=ind+1;
    end
end
m_train=m_train(v_usersToKeep,:);
m_test=m_test(v_usersToKeep,:);
c_userInfo=c_userInfo(v_usersToKeep,:);
c_movieInfo=c_movieInfo(v_usersToKeep,:);
v_diffentOc=unique(c_userInfo(:,4));
v_diffentZipCod=unique(c_userInfo(:,5));
v_diffentSex=unique(c_userInfo(:,3));
v_ageGroups=[0,12;13,15;16,17;18,20;21,24;25,28;29,33;34,38;39,42;43,46;47,50;51,55;56,60;61,65;66,200];
v_randMixingVector=[0.5,0.7,0.9];
m_clustUserInfoOccup=zeros(size(c_userInfo,1),size(v_diffentOc,1));
%m_clustUser contains all the cluster indicators
m_clustUser=zeros(size(c_userInfo,1),size(v_diffentOc,1)+size(v_diffentSex,1)+size(v_ageGroups,1));
for k=1 :size(v_diffentOc,1)
    m_clustUser(:,k)=strcmp(c_userInfo(:,4),v_diffentOc(k))*v_randMixingVector(1);
end
m_clustUserInfoSex=zeros(size(c_userInfo,1),size(v_diffentSex,1));

for k=1 :size(v_diffentSex,1)
    m_clustUser(:,k+size(v_diffentOc,1))=strcmp(c_userInfo(:,3),v_diffentSex(k))*v_randMixingVector(2);
end
%%Cluster ages
%based on your intuition 
%based on the data so that it is more fine grained
m_clustUserInfoAge=zeros(size(c_userInfo,1),size(v_ageGroups,1));
v_ages=cell2mat(c_userInfo(:,2));
for k=1 :size(v_ageGroups,1)
    m_clustUser(:,k+size(v_diffentSex,1)+size(v_diffentOc,1))=((v_ageGroups(k,1)<=v_ages)&(v_ages<=v_ageGroups(k,2)))*v_randMixingVector(3);
end



m_clustUserInfoAge=m_clustUser(:,end-size(v_ageGroups,1)+1:end);
m_clustUserInfoOccup=m_clustUser(:,1:size(v_diffentOc,1));
m_clustUserInfoSex=m_clustUser(:,size(v_diffentOc,1)+1:end-size(v_ageGroups,1));
%add random generated row.... as a last row
%m_clustUser=m_
end