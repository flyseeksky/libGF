
classdef ReadMovieLensDataset < ReadDataset
	
	
	properties(Constant)
		ch_folderName = './libGF/datasets/MovieLensDataset/ml-100k/';
		
	end
	
	
	
	methods(Static)
		
		function [ m_test, v_range ] = getTestTables			
			% Obtain 5 disjoint sets  with 20.000 ratings  to reproduce the
			% simulations in [narang2013structured]
			%
			% m_test        : 943 x 1682 x 5 array. m_test(:,:,i) contains
			%               the  i-th test set. It has one row per user and
			%               one column per movie. 
		    %               
			%               The i-th slab is obtained from ui.data
			%
			% v_range       : 2x1 vector contains the min and max rating 
			%               available in the dataset
			%
				
			s_SetNum = 5;
			s_UsersNum = 943; 
			s_MoviesNum = 1682;			
			m_test=zeros(s_UsersNum,s_MoviesNum,s_SetNum);
			v_range=zeros(2,1);
			
			for s_setInd = 1:s_SetNum
				
				ch_file = sprintf('%s%d%s',[ReadMovieLensDataset.ch_folderName 'u'],s_setInd,'_test.txt');
				
				D = readtable(ch_file,'Delimiter','\t','ReadVariableNames',false);
				TestData=table2array(D); % columns contain user id | item id | rating | timestamp.
				%                          The time stamps are unix seconds since 1/1/1970 UTC
				
				for icount=1:size(TestData,1)
					m_test(TestData(icount,1),TestData(icount,2),s_setInd)=TestData(icount,3);
				end
			
			end		
			%m_test=(m_test-1)/4;
			m_test(m_test==0)=NaN;
			v_range(1)=min(min(min(m_test)));
			v_range(2)=max(max(max(m_test)));
			% NaN and range [0,1]
			
			%save([ReadMovieLensDataset.ch_folderName 'm_test.mat'],'m_test','v_range');
		end
		
	end
end
