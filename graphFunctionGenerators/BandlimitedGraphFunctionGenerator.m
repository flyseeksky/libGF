classdef BandlimitedGraphFunctionGenerator  < GraphFunctionGenerator
	
	
	properties % Required by superclass Parameter
		c_parsToPrint    = {'ch_name','s_bandwidth'};
		c_stringToPrint  = {'',    'B'};
		c_patternToPrint = {'%s%s signal','%s = %d'};
	end 
	
	properties
		ch_name = 'Bandlimited';
		s_bandwidth    % integer
		
		b_sortedSpectrum = 0;  % if 1, the entries of the Fourier transform 
		% are generated and then sorted
        
        b_generateSameFunction = 0; % generate the same function if set to 1
	
		
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
			
            if obj.b_generateSameFunction
                rng(1);
            end
			
			m_B = obj.basis;
            %freq = randn(obj.s_bandwidth,s_numberOfRealizations);
			if obj.b_sortedSpectrum
                %[~, ind] = sort(abs(freq), 'descend');
				%M_graphFunction = sqrt(size(m_B,1)/obj.s_bandwidth) * m_B*freq(ind);
                M_graphFunction = sqrt(size(m_B,1)/obj.s_bandwidth) * ...
					m_B*sort(randn(obj.s_bandwidth,s_numberOfRealizations),1,'descend');
			else
				M_graphFunction = sqrt(size(m_B,1)/obj.s_bandwidth) * ...
					m_B*randn(obj.s_bandwidth,s_numberOfRealizations);
			end
				
				%             N = obj.graph.getNumberOfVertices();
				%             atilde = (1:N)';
				%             alpha = exp( - atilde / 50 );
				%             V = obj.graph.getLaplacianEigenvectors();
				%             M_graphFunction = 10*V*alpha * ones(s_numberOfRealizations,1);
				
			
		end
		
		function m_basis = basis(obj,s_otherBandwidth)
			%  M_BASIS            N x S_OTHERBANDWIDTH matrix containing the
			%                     first OBJ.s_bandwidth eigenvectors of the
			%                     Laplacian
			%  S_OTHERBANDWIDTH   optional parameter to specify the number of 
			%                     desired columns. (Default: =
			%                     obj.s_bandwidth)
			
			if nargin<2 % default option
				s_otherBandwidth = obj.s_bandwidth;
			end
			m_V = obj.graph.getLaplacianEigenvectors();
			m_basis = m_V(:,1:s_otherBandwidth);
			
		end
		
		
	end
	
end

