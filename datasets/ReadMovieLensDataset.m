
classdef ReadMovieLensDataset
	
	
	methods
		
		
	end
	
	
	
	methods(Static)
		
		function [ m_test ] = getTestTables			
			% Obtain 5 disjoint sets  with 20.000 ratings  to reproduce the
			% simulations in [narang2013structured]
			%
			% m_test        : 943 x 1682 x 5 array. m_test(:,:,i) contains
			%               the  i-th test set. It has one row per user and
			%               one column per movie. 
		    %               
			%               The i-th slab is obtained from ui.data
			
					
			ch_folderName = './libGF/datasets/MovieLensDataset/ml-100k/';
			
			s_SetNum = 5;
			s_UsersNum = 943; 
			s_MoviesNum = 1682;			
			m_test=zeros(s_UsersNum,s_MoviesNum,s_SetNum);
			
			
			for s_setInd = 1:s_SetNum
				
				ch_file = sprintf('%s%d%s',[ch_folderName 'u'],s_setInd,'_test.txt');
				
				D = readtable(ch_file,'Delimiter','\t','ReadVariableNames',false);
				TestData=table2array(D); % columns contain user id | item id | rating | timestamp.
				%                          The time stamps are unix seconds since 1/1/1970 UTC
				
				for icount=1:size(TestData,1)
					m_test(TestData(icount,1),TestData(icount,2),s_setInd)=TestData(icount,3);
				end
			
			end			
			
		end
		
	end
end
