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
				
		m_estimate = estimate(obj,m_samples,m_positions);			
		%
		% Input:
		% M_SAMPLES                 S x S_NUMBEROFREALIZATIONS  matrix with
		%                           samples of the graph function in
		%                           M_GRAPHFUNCTION 
		% M_POSITIONS               S x S_NUMBEROFREALIZATIONS matrix
		%                           containing the indices of the vertices
		%                           where the samples were taken
		%
		% Output:                   
		% M_ESTIMATE                N x S_NUMBEROFREALIZATIONS matrix. N is
		%                           the number of nodes and each column
		%                           contains the estimate of the graph
		%                           function
		% 
		
		
		
	end
	
end

