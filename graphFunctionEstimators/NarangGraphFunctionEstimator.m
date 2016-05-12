classdef NarangGraphFunctionEstimator < GraphFunctionEstimator	
	% Estimators from [narang2013localized]
	properties(Constant)
	end
	
	properties % Required by superclass Parameter
		c_parsToPrint    = {'ch_name','ch_type'};
		c_stringToPrint  = {'',''};
		c_patternToPrint = {'%s%s','%s%s'};
	end
	
	properties
		ch_name = 'Narang et al. 2013';
		
		ch_type = 'RBM'; % it can be
		%      'RBM':    regularization-based method
		%      'LSR':    least-squares reconstruction
				
		s_theta = 1; % parameter for Bilateral Link-Weight Adjustment
		
	end
	
	properties(Access = private) % precomputed properties for efficiency
		v_laplacianEigenvalues_precomp % sorted in increasing order
		m_laplacianEigenvectors_precomp
	end
	
	methods
		
		function obj = NarangGraphFunctionEstimator(varargin)
			obj@GraphFunctionEstimator(varargin{:});
		end
		
		
		
		function estimate = estimate(obj,m_samples,sideInfo)
			%
			% Input:
			% M_SAMPLES                 S x S_NUMBEROFREALIZATIONS  matrix with
			%                           samples of the graph function in
			%                           M_GRAPHFUNCTION
			% sideInfo
			%         1 x S_NUMBEROFREALIZATIONS vector of structs with fields
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
			% estimate
			%         1 x S_NUMBEROFREALIZATIONS vector of structs with
			%         fields
			%         estimate(i).v_wantedSamples: W x 1 vector containing the
			%                            estimated signal at the entries
			%                            indicated by
			%                            sideInfo(i).v_wantedEntries
			%
			%
								
			% input check
			if size(m_samples,2)>1
				error('not implemented')				
			end
			assert( isa(sideInfo,'struct') );
			assert( isfield(sideInfo,'v_sampledEntries') );
			assert( isrow(sideInfo.v_sampledEntries') );
			assert( isfield(sideInfo,'graph') );
			s_numberOfVertices = sideInfo.getNumberOfVertices(); % number of vertices
			if  ~isfield(sideInfo,'v_wantedEntries') 
				sideInfo.v_wantedEntries = (1:s_numberOfVertices)';
			end
			assert( isrow(sideInfo.v_wantedEntries') );
			
			% compute subgraph with sampled and wanted vertices
			m_graphAdjacency = sideInfo.graph.m_adjacency;
			v_graphIndices = [sideInfo.v_sampledEntries;sideInfo.v_wantedEntries];
			m_subgraphAdjacency = m_graphAdjacency( v_graphIndices , v_graphIndices );
			subgraph = Graph('m_adjacency',m_subgraphAdjacency);			
			
			% nearest neighbors graph
			subgraph = subgraph.nearestNeighborsSubgraph(30);
			
			% update weights			
			subgraph = obj.bilateralLinkWeightAdjustment( subgraph , 1:length(sideInfo.v_sampledEntries) , m_samples , obj.s_theta ); 
			
			% estimate subsignal
			estimate.v_wantedSamples = obj.estimateSubsignal(subgraph, 1:length(sideInfo.v_sampledEntries) , m_samples , 1+length(sideInfo.v_sampledEntries):length(v_graphIndices) );
			
			
			%-------------------------------------------------------------
			% 
% 			s_numberOfVertices = size(obj.m_laplacian,1);
% 			s_numberOfRealizations = size(m_samples,2);
% 			m_estimate = zeros(s_numberOfVertices,s_numberOfRealizations);
% 			
% 			
% 			for iRealization = 1:s_numberOfRealizations
% 				
% 				if obj.s_bandwidth == -1
% 					s_bw = obj.computeCutoffFrequency(m_positions(:,iRealization));
% 					m_eigenvecs = obj.m_laplacianEigenvectors_precomp(:,1:s_bw);
% 				end
% 				
% 				m_PhiB = m_eigenvecs( m_positions(:,iRealization) , : );
% 				v_alphas = m_PhiB\m_samples(:,iRealization);
% 				m_estimate(:,iRealization) = m_eigenvecs*v_alphas;
% 			end
			
			
		end
			
		
		function s_bw = computeCutoffFrequency(obj, m_positions)
			
			
			L = obj.m_laplacian;
			L_power = L^(2*obj.proxy_order);
			
			% svds(L_power(m_positions,m_positions),1,0) does not always work, even
			% after setting a high tolerance
			v_svals = svd(L_power(m_positions,m_positions));
			min_sval = min(v_svals);
			omega_S = (  min_sval  )^(1/(2*obj.proxy_order));
			
			
			s_bw = sum(obj.v_laplacianEigenvalues_precomp <= omega_S);
			
		end
		
		
		
		function subsignal = estimateSubsignal(obj,graph, v_sampledEntries , m_samples , v_wantedEntries )
			
			switch obj.ch_type
				case 'RBM'
					
				case 'LSR'
					
			end
			
		end
		
		
	end
	methods(Static)
		
		function graph = bilateralLinkWeightAdjustment( graph , v_entries , v_samples , s_theta )
			% v_samples(i) is the value of the graph function on vertex 
			% v_entries(i)
			%
			% output graph is a modified graph using [narang2013structured, eq
			% (10)]
			%
			W = graph.m_adjacency;
            for row = 1 : size(W,1)
                for col = 1 : size(W,1)
                    if W(row, col) > 0
                        W(row, col) = W(row,col) * exp( - (v_samples(row) - ...
                            v_samples(col) ) / s_theta );
                    end
                end
            end
            graph = Graph('m_adjacency', W);				
        end
	end

end
