classdef RecommenderSystemsSimulator < Simulator
	
	properties
	end
	
	methods(Static)
		
		function mse = simulateDataset( v_CVSets , graphConstructor, v_estimator )
			% v_CVSets:    L x 1 vector of structs with fields
			%     -v_CVSets(i).m_training
			%     -v_CVSets(i).m_validation
			%    Both these two fields are moviesNum x usersNum matrices,
			%    where the (i,j) element contains the rating of user j to
			%    movie i. Misses (unknown entries) are marked with NaN.
			%
			% graphConstructor: function that, given a matrix like
			%    v_CVSets(i).m_training, it returns an object of class graph.
			%
			% v_estimator:   E x 1 vector of objects of class GraphFunctionEstimator
			%
			% mse is computed per estimator v_estimator(k) as follows:
			%    1. v_mse(i) is computed for the i-th set, i.e., v_CVSets(i),
			%    for all i. To do this
			%          1- a graph is constructed using graphConstructor
			%          applied to v_CVSets(i).m_training
			%          2- estimator is invoked per column to estimate the
			%          non-empty (different from NaN) entries in
			%          v_CVSets(i).m_validation 
			%          3- v_mse(i) is computed over the non-empty entries of
			%          v_CVSets(i).m_validation 
			%    2. mse is the result of averaging v_mse
			%             
			%         

			

			s_movieNum = size(v_CVSets(1).m_training,1);			
			s_userNum = size(v_CVSets(1).m_training,2);			
			s_foldCVNum = length(v_CVSets);
			s_estimatorNum = length(v_estimator);
			
			m_mse = NaN(s_estimatorNum,s_foldCVNum);
			v_rowIndices = (1:s_movieNum)';
			


			for s_foldCVInd = 1:s_foldCVNum
				s_foldCVInd
load('dummy.mat')				
%				graph_now = graphConstructor(v_CVSets(s_foldCVInd).m_training);
%save('dummy.mat')
                
                % update estimators
				for s_estimatorInd = s_estimatorNum:-1:1
					v_estimator(s_estimatorInd) = v_estimator(s_estimatorInd).prepareForGraph(graph_now);					
				end

				for s_userInd = s_userNum:-1:1
					% estimation of missing entries column by column
					s_userInd
					v_training = v_CVSets(s_foldCVInd).m_training(:,s_userInd);
					v_trainingEntries = v_rowIndices(~isnan(v_training));
					v_trainingSamples = v_training(v_trainingEntries);
					
					v_validation = v_CVSets(s_foldCVInd).m_validation(:,s_userInd);
					v_validationEntries = v_rowIndices(~isnan(v_validation));
					v_validationSamples = v_validation(v_validationEntries);
					
					% estimation
					sideInfo.v_sampledEntries = v_trainingEntries;
					sideInfo.v_wantedEntries = v_validationEntries;					

					for s_estimatorInd = s_estimatorNum:-1:1
						v_signalEstimate = v_estimator(s_estimatorInd).estimate(v_trainingSamples,sideInfo);
						
						% error computation
						v_squaredError(s_estimatorInd,s_userInd) = norm( v_validationSamples - v_signalEstimate.v_wantedSamples )^2;
						s_estimatedEntriesNum(s_estimatorInd,s_userInd) = length( v_validationEntries );
					end
					
				end
				m_mse(:,s_foldCVInd) = sum( v_squaredError , 2 )./sum( s_estimatedEntriesNum , 2 );
				%end
			end
			m_mse
			mse = mean(m_mse,2);
			
			
		end
		
		
	end
	
end

