classdef GraphLearningSmoothSignalGraphGenerator < GraphGenerator
    properties % required by parent classes
        c_parsToPrint  = {};
        c_stringToPrint  = {};
        c_patternToPrint = {};
    end
    
    properties(Constant)
        ch_name = 'Graph-Learning-For-Smooth-Signal';
    end
    
    properties
        m_observed; %contains the observed signals
        s_niter;    %maximum number of iterations
        s_alpha;
        s_beta; %positive regularization parameters0000......
        
        m_missingValuesIndicator;   %is a matrix with one if the corresponding entry in m_observed is
                                    %observed and 0 otherwise. If [] is used then assumed that all values observed
    end
    
    
    methods
        function obj = GraphLearningSmoothSignalGraphGenerator(varargin)
            % Constructor
            obj@GraphGenerator(varargin{:});
            
        end
        
        function graph = realization(obj)
            L = obj.learnLaplacian(obj.m_observed, obj.s_alpha, obj.s_beta);              
            m_adjacency=Graph.createAdjacencyFromLaplacian(L);
            graph = Graph('m_adjacency',m_adjacency);
        end
    end
    
    
    methods (Static)
        function m_laplacian = learnLaplacian(m_observed, s_alpha, s_beta)
            % learning laplacian matrix L by alternating minimization
            % this is the first step, i.e. fix Y minimize w.r.t L
            % 
            N = size(m_observed,1); % number of nodes
            Y = m_observed;         % observed signals, in columns
            M = GraphLearningSmoothSignalGraphGenerator.getMdup(N);
            m_A = GraphLearningSmoothSignalGraphGenerator.getA(N);
            m_B = GraphLearningSmoothSignalGraphGenerator.getB(N);
            s_length = N*(N+1)/2;
            
            cvx_begin quiet
                variable vech_L(s_length) 
                minimize( s_alpha * vec(Y*Y')' * M * vech_L + ...
                    s_beta * vech_L' * M' * M * vech_L )
                subject to
                    m_A * vech_L == N;
                    m_B * vech_L <= 0;
            cvx_end
            
            m_laplacian = GraphLearningSmoothSignalGraphGenerator.ivech(vech_L);
        end
        
        function m_A = getA(N)
            m_A = zeros(N, N*(N+1)/2);
            v_diagIndex = GraphLearningSmoothSignalGraphGenerator.getDiagIndexInVech(N);
            for row = 1 : N
                col = v_diagIndex(row);
                m_A(row, col) = 1;
            end
        end
        
        function m_B = getB(N)
            % get matrix B
            % pick all the off-diagonal entries from vech(L)
            m_B = zeros(N*(N-1)/2, N*(N+1)/2);
            diagIndex = GraphLearningSmoothSignalGraphGenerator.getDiagIndexInVech(N);
            row = 1;
            for idx = 1 : N*(N+1)/2
                if any(idx == diagIndex)
                    continue;
                else
                    m_B(row, idx) = 1;
                    row = row + 1;
                end
            end
        end
        
        function v_diagIndex = getDiagIndexInVech(N)
            % for a graph of given size N
            % find the index of diagonal elements in vector vech(L)
            v_diagIndex = 1 + [0 cumsum(N:-1:2)];
        end
        
        function M = getMdup(s_laplacianSize)
            % generate duplication matrix M such that
            % M * vech(L) = vec(L)
            N = s_laplacianSize;
            vSize = N^2;
            vhSize = N*(N+1)/2;
            M = zeros(vSize,vhSize);
            
            % convert index of vec(L) into coordinates
            % then convert coordinates into index of vech(L)
            colIndex = [0 cumsum(N:-1:1)];
            for idx = 1:vSize
                coord = [mod(idx-1, N)+1, floor((idx-1)/N)+1];
                row = max(coord);
                col = min(coord);
                idxh = colIndex(col) + row - col + 1;
                M(idx, idxh) = 1;
            end
        end
        
        function m_matrixData=ivech(v_stackedData)
            % converts the vectorized form of lower triangular part of a matrix back
            v_stackedData = v_stackedData(:);
            len = length(v_stackedData);
            N = (-1+sqrt(1+8*len))/2;
            
            if floor(N)~=N
                error(['The number of elemeents in STACKED_DATA must be conformable to' ...
                    'the inverse vech operation.'])
            end
            
            m_matrixData=zeros(N);
            
            pl=tril(true(N));
            m_matrixData(pl)=v_stackedData;
            diag_matrixData=diag(diag(m_matrixData));
            m_matrixData=m_matrixData+m_matrixData'-diag_matrixData;
        end
        
        function v_vech = vech(m_X)
            % contains the vectorized form of lower triangular part of m_X
            % m_Y=m_X';
            v_vech = m_X(tril(true(size(m_X))));
        end
    end
end