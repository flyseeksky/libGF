function [ m_averagedTrainData] = handleMLdat(m_train,m_test,m_userInfo,m_movieInfo)
%As signals all the movies ?
%Or for less missing data cluster the movies and take averaging rating of a
%user for each category  This way use as training data only the averaged
%signals tha result from the caregories
catIndicatorVectors=cell2mat(m_movieInfo(:,7:end));
movOfCat=zeros(size(m_movieInfo,1),size(catIndicatorVectors,2));
m_averagedTrainData=zeros(size(m_train,1),size(catIndicatorVectors,2));
for i=1:size(catIndicatorVectors,2)
    a=cell2mat(m_movieInfo(catIndicatorVectors(:,i)==1,1));
    movOfCat(1:size(a,1),i)=a;
    b=sum(m_train(:,a),2);
    nonzeroRatings=sum(m_train(:,a)~=0,2);
    m_averagedTrainData(:,i)=b./nonzeroRatings;
    
end
m_averagedTrainData(isnan(m_averagedTrainData))=0;
%the columns of AveragedTrainData contain the average rating of each user
%for every kind o movie
end

