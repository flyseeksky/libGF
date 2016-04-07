classdef Graph
	% Contains a graph and allows to perform several operations
	%
	
	properties
		m_adjacency % adjacency matrix
	end
	
	methods
				
		function obj = Graph(varargin)
		% Constructor	
			if nargin > 0
				obj =  assignParametersByName(obj,varargin{:});
			end	
		end
		
		
		function m_L = getLaplacian(obj)						
			v_degrees = sum(obj.m_adjacency,2);
			m_L = diag(v_degrees) - obj.m_adjacency;
		end
		
		function m_V = getLaplacianEigenvectors(obj)									
			% m_V     Matrix whose columns are the eigenvectors of the
			%         Laplacian sorted in ascending order of eigenvalue
			
			m_L = obj.getLaplacian;
			[m_V,~] = eig(m_L);
		end
		
		
		function s_n = getNumberOfVertices(obj)
			s_n = size(obj.m_adjacency,1);			
        end
        

		
		function c_components = getComponents(obj)
			% (To be written)
			%
			% COMPONENTS is a cell array of C vectors, where C is the
			% number of components of the graph OBJ. COMPONENTS{c} is a
			% column vector containing the indices of the vertices in each
			% component. 
			%
			% The algorithm used...
            m_sparseAdjacency=sparse(obj.m_adjacency);
            [s_numberOfComponents,v_componentIndicator]=graphconncomp(m_sparseAdjacency,'DIRECTED',false,'WEAK',true);
            %v_componentIndicator is a vector that for each vertice has the
            %corresponding component number 
            %in s_numberOfComponents is the number of components
			c_components = {};
            for s_k=1:s_numberOfComponents
                c_components(:,s_k)={(find(v_componentIndicator(:)==s_k))};
            end
        end
		
		function C = getClusters(obj,s_numberOfClusters,s_Type)
            % Executes the spectral clustering algorithm defined by
            %   Type on the adjacency matrix W and returns the k cluster
            %   indicator vectors as columns in C.
            %   If L and U are also called, the (normalized) Laplacian and
            %   eigenvectors will also be returned.
            %
            % Input:
            % S_NUMBEROFCLUSTERS   Number of clusters
            %   'Type' - Defines the type of spectral clustering algorithm
            %            that should be used. Choices are:
            %      1 - Unnormalized
            %      2 - Normalized according to Shi and Malik (2000)
            %      3 - Normalized according to Jordan and Weiss (2002)
            %
            %
			% Output: 
			% C          
            %
            %   References:
            %   - Ulrike von Luxburg, "A Tutorial on Spectral Clustering",
            %     Statistics and Computing 17 (4), 2007
            %
            
            W = obj.m_adjacency;
            k = s_numberOfClusters;
            
            % calculate degree matrix
            degs = sum(W, 2);
            D    = sparse(1:size(W, 1), 1:size(W, 2), degs);
            
            % compute unnormalized Laplacian
            L = D - W;
            
            % compute normalized Laplacian if needed
            switch s_Type
                case 2
                    % avoid dividing by zero
                    degs(degs == 0) = eps;
                    % calculate inverse of D
                    D = spdiags(1./degs, 0, size(D, 1), size(D, 2));
                    
                    % calculate normalized Laplacian
                    L = D * L;
                case 3
                    % avoid dividing by zero
                    degs(degs == 0) = eps;
                    % calculate D^(-1/2)
                    D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
                    
                    % calculate normalized Laplacian
                    L = D * L * D;
            end
            
            % compute the eigenvectors corresponding to the k smallest
            % eigenvalues
            diff   = eps;
            [U, ~] = eigs(L, k, diff);
            
            % in case of the Jordan-Weiss algorithm, we need to normalize
            % the eigenvectors row-wise
            if s_Type == 3
                U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
            end
            
            % now use the k-means algorithm to cluster U row-wise
            % C will be a n-by-1 matrix containing the cluster number for
            % each data point
            C = kmeans(U, k, 'start', 'cluster', ...
                'EmptyAction', 'singleton');
            
            % now convert C to a n-by-k matrix containing the k indicator
            % vectors as columns
            C = sparse(1:size(D, 1), C, 1);
            
        end
        
				
	end
	
	methods(Static)
% 		function G = constructViaFactorAnalysis( m_functionValues , alpha , beta )
% 			% (To be written)
% 			%
% 		    % Construct graph from signal values using [Dong et al. 2015]
% 			%
% 			% Input:
% 			% M_FUNCTIONVALUES     N x M Matrix where N is the number of
% 			%                      nodes and M is the number of
% 			%                      observations of a function on a graph
% 			% ALPHA, BETA          Regularization parameters of the
% 			%                      algorithm
% 			%
% 			% Output:
% 			% G                    Graph of N nodes
% 			%
% 			
% 			m_adjacency_est = [];
% 			G = Graph('m_adjacency',m_adjacency_est);
% 			
% 		end
		function m_adjacency= createAdjacencyFromLaplacian(m_laplacian)
            m_adjacency=eye(size(m_laplacian)).*diag(diag(m_laplacian))-m_laplacian;

        end
		
	end
	
	
end

