classdef BandlimitedGraphFunctionEstimator < GraphFunctionEstimator
	% 
	properties(Constant)
	end
	
	properties % Required by superclass Parameter
		c_parsToPrint    = {};
		c_stringToPrint  = {};
		c_patternToPrint = {};
	end
	
	properties
		ch_name = 'BANDLIMITED';
		m_laplacianEigenvectors;   % N x s_bandwidth matrix with the first
		% s_bandwidth eigenvectors of the
		% Laplacian
		
	end
	
	methods
		
		function obj = BandlimitedGraphFunctionEstimator(varargin)
			obj@GraphFunctionEstimator(varargin{:});
		end
		
	end
	
	methods
		
		function m_estimate = estimate(obj,m_samples,m_positions)
			%
			% Input:
			% M_SAMPLES                 S x S_NUMBEROFREALIZATIONS  matrix with
			%                           samples of the graph function in
			%                           M_GRAPHFUNCTION
			% M_POSITIONS               S x S_NUMBEROFREALIZATIONS matrix
			%                           containing the indices of the vertices
			%                           where the samples were taken
			%
			% Output:                   N x S_NUMBEROFREALIZATIONS matrix. N is
			%                           the number of nodes and each column
			%                           contains the estimate of the graph
			%                           function
			%
			
			s_numberOfVertices = size(obj.m_laplacianEigenvectors,1);
			s_numberOfRealizations = size(m_samples,2);
			
			m_estimate = zeros(s_numberOfVertices,s_numberOfRealizations);
			for realizationCounter = 1:s_numberOfRealizations
				m_PhiB = obj.m_laplacianEigenvectors( m_positions(:,realizationCounter) , : );
				v_alphas = m_PhiB\m_samples(:,realizationCounter);
				m_estimate(:,realizationCounter) = obj.m_laplacianEigenvectors*v_alphas;
			end
			
			
		end
		
	end

end
