classdef BandlimitedGraphFunctionGenerator  < GraphFunctionGenerator
	
	
	properties % Required by superclass Parameter
		c_parsToPrint    = {'name','s_bandwidth'};
		c_stringToPrint  = {'',    'BANDWIDTH'};
		c_patternToPrint = {'%s%s','%s = %d'};
	end 
	
	properties
		name = 'Bandlimited';
		s_bandwidth    % integer		
	end
	
	methods
		
		function obj = BandlimitedGraphFunctionGenerator(varargin)
			% constructor
			obj@GraphFunctionGenerator(varargin{:});
		end
		
		
		function M_graphFunction = realization(obj,s_numberOfRealizations)
			% M_GRAPHFUNCTION   N x S_NUMBEROFREALIZATIONS matrix where N is
			%                   the number of vertices. Each column is a
			%                   signal whose graph fourier transform is 
			%                   i.i.d. standard Gaussian distributed for
			%                   the first OBJ.s_bandwidth entries and zero
			%                   for the remaining ones
			
			assert(~isempty(obj.graph));
			assert(~isempty(obj.s_bandwidth));
			
			if nargin < 2
				s_numberOfRealizations = 1;
			end
			
			
			m_B = obj.basis;
			M_graphFunction = sqrt(size(m_B,1)/obj.s_bandwidth)*m_B*randn(obj.s_bandwidth,s_numberOfRealizations);
			
			
		end
		
		function m_basis = basis(obj)
			%  M_BASIS       N x OBJ.s_bandwidth matrix containing the
			%                first OBJ.s_bandwidth eigenvectors of the
			%                Laplacian
			
			m_V = obj.graph.getLaplacianEigenvectors();
			m_basis = m_V(:,1:obj.s_bandwidth);
			
		end
		
		
	end
	
end

