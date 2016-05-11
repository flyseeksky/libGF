classdef GraphFunctionEstimator < Parameter
	% This class is cool
	properties(Constant)
	end
	
	properties
	end
		
	methods
		
		function obj = GraphFunctionEstimator(varargin)
			obj@Parameter(varargin{:});
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

