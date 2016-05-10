classdef UnknownBandwidthBandlimitedGraphFunctionEstimator < GraphFunctionEstimator
	% This is a comment
	properties(Constant)
	end
	
	properties % Required by superclass Parameter
		c_parsToPrint    = {'ch_name','s_bandwidth'};
		c_stringToPrint  = {'','ASS. BW'};
		c_patternToPrint = {'%s%s','%s = %d'};
	end
	
	properties
		ch_name = 'BUB';
		m_laplacian;  % N x N Laplacian matrix
			
	end
	
	methods
		
		function obj =  UnknownBandwidthBandlimitedGraphFunctionEstimator(varargin)
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
			
			if ~isempty(obj.s_bandwidth)
				if( obj.s_bandwidth > size(obj.m_laplacianEigenvectors,2) )
					error('s_bandwidth cannot be greater than the number of columns provided in m_laplacianEigenvectors');
				end
				obj.m_laplacianEigenvectors = obj.m_laplacianEigenvectors(:,1:obj.s_bandwidth);
			end
			
			s_numberOfVertices = size(obj.m_laplacianEigenvectors,1);
			s_numberOfRealizations = size(m_samples,2);
						
			m_estimate = zeros(s_numberOfVertices,s_numberOfRealizations);
			for iRealization = 1:s_numberOfRealizations
				m_PhiB = obj.m_laplacianEigenvectors( m_positions(:,iRealization) , : );
				v_alphas = m_PhiB\m_samples(:,iRealization);
				m_estimate(:,iRealization) = obj.m_laplacianEigenvectors*v_alphas;
			end
			
			
		end
		
	end

end
