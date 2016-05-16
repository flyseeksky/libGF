classdef GraphFunctionEstimator < Parameter
	% This class is cool
	properties(Constant)
	end
	
	properties
		s_regularizationParameter
		s_numFoldValidation
	end
		
	methods
		
		function obj = GraphFunctionEstimator(varargin)
			obj@Parameter(varargin{:});
		end
		
		
		function [s_optMu,s_ind] = crossValidation(obj,v_samples,v_positions,v_mu)
			% Input:
			% V_SAMPLES                 S x S_NUMBEROFREALIZATIONS  matrix with
			%                           samples of the graph function in
			%                           M_GRAPHFUNCTION
			% V_POSITIONS               an S x S_NUMBEROFREALIZATIONS
			%                           matrix containing the indices of
			%                           the vertices where the samples were
			%                           taken  
			
			assert(size(v_samples,2)==1,'not implemented');
			
			m_mse = zeros(1,length(v_mu));
			for muInd = 1:length(v_mu)
				
				% partition v_positions in obj.s_numFoldValidation subsets
				% m_cvPositions =   % N0 x obj.s_numFoldValidation matrix
				
				for valInd = 1:obj.s_numFoldValidation
					
					% Create test and validation set
					% v_samples_training =
					% v_positions_training =
					% v_samples_validation =
					% v_positions_validation =    % indices of the "unobserved"
					%                       % vertices in crossvalidation
					
					
					% Estimate
					obj.s_regularizationParameter = v_mu(muInd);
					v_signal_est = obj.estimate(obj,m_samples_training,v_positions_training);
					
					% Measure MSE
					
					m_mse(muInd,valInd) = norm( v_samples(v_positions_validation) - v_signal_est(v_positions_validation) )
				end
				
				
			end
			v_mse = mean(m_mse,2);
			[~,s_ind] = min(v_mse);
			s_optMu = v_mu(s_ind);
		end
		
	end
	
	methods(Abstract)
				
		estimate = estimate(obj,m_samples,sideInfo);			
		%
		% Input:
		% M_SAMPLES                 S x S_NUMBEROFREALIZATIONS  matrix with
		%                           samples of the graph function in
		%                           M_GRAPHFUNCTION 
		% sideInfo                  It can be either:
		%      a) an S x S_NUMBEROFREALIZATIONS matrix containing the
		%      indices of the vertices where the samples were taken
		%      b) a 1 x S_NUMBEROFREALIZATIONS vector of structs with fields
		%         sideInfo(i).v_sampledEntries:  S x 1
		%                           vector where each column contains the
		%                           indices of the sampled vertices
		%         sideInfo(i).v_wantedEntries:  W x 1
		%                           vector where each column contains the
		%                           indices of the desired vertices. If not
		%                           defined, it is assumed that this field
		%                           is 1:N.
		%         sideInfo(i).graph: graph over which the signal has been
		%                           sampled
		%
		% Output:                   
		% estimate                  It can be either:
		%      a) an N x S_NUMBEROFREALIZATIONS matrix. N is the number of
		%      nodes and each column contains the estimate of the graph
		%      function 
		%      b) a 1 x S_NUMBEROFREALIZATIONS vector of structs with
		%      fields
		%         estimate(i).v_wantedSamples: W x 1 vector containing the
		%                            estimated signal at the entries
		%                            indicated by
		%                            sideInfo(i).v_wantedEntries 
		%          
		% 
		
		
		
	end
	
end

