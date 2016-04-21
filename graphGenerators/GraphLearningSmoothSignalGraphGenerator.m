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
        
        function [graph,m_estimated] = realization(obj)
            %initialization
            m_missingValuesIndicator=obj.m_missingValuesIndicator;
            m_observed=obj.m_observed;
            s_alpha=obj.s_alpha;
            s_beta=obj.s_beta;
            s_niter=obj.s_niter;
            
            if(isempty(m_missingValuesIndicator)) %#ok<*PROP>
                m_missingValuesIndicator=ones(size(m_observed));
            end
            m_estimated=m_observed;
            %m_dublication converts the vectorized form of lower triangular
            %part of L to the vectorized form of the full matrix
            m_dubplication=GraphLearningSmoothSignalGraphGenerator.dup(rand(size(m_estimated,1)));
            
            
            %alternating optimization
            v_obj=zeros(s_niter,1);
            for t=1:s_niter
                m_laplacian=GraphLearningSmoothSignalGraphGenerator.graphLaplUpd(s_alpha,s_beta,m_estimated,m_dubplication);
                m_estimated=GraphLearningSmoothSignalGraphGenerator.signalUpd(m_observed,m_laplacian,s_alpha,m_missingValuesIndicator);
                %contains the objective values of the function
                v_obj(t)=(norm(m_missingValuesIndicator.*(m_observed-m_estimated),'fro')^2)+s_alpha*trace((m_estimated)'*m_laplacian*(m_estimated))+s_beta*norm(m_laplacian,'fro')^2;
                if t>1&&(0<(-v_obj(t)+v_obj(t-1)))&&((-v_obj(t)+v_obj(t-1))<10^-6)
                    break;
                end
            end
            % figure(1);
            % semilogy(1:t,v_obj);
            %keep only meaningful values
            %Prune insignificant edges
            s_alpha=(m_laplacian>-10^-4);
            s_beta=m_laplacian>10^-4;
            s_alpha=~s_alpha;
            c=s_alpha|s_beta;
            m_laplacian(~c)=0;
                          
            m_adjacency=Graph.createAdjacencyFromLaplacian(m_laplacian);
            graph = Graph('m_adjacency',m_adjacency);
        end
        
        function m_L = test(obj)
            m_L = obj.learnLaplacian(obj.m_observed, obj.s_alpha, obj.s_beta);
        end
    end
    
    
    methods (Static)
        %f
        function m_laplacian=graphLaplUpd(s_alpha,s_beta,m_observed,m_duplication)
            % arg min a*vec(Y*Y')'*Mdup*vech(L)+b*vech(L)'*Mdup'*Mdup*vech(L)
            % wrt vech(L)
            % st A*vech(L)=0;
            %   B*vech(L)<=0;
            k=size(m_observed,1);
            % Is used for the equality constraints of the opt problem
            m_A=GraphLearningSmoothSignalGraphGenerator.getLinEqualities(rand(size(m_observed,1)));
            m_B=GraphLearningSmoothSignalGraphGenerator.getLinInequalities(rand(size(m_observed,1)));
            s_n=size(GraphLearningSmoothSignalGraphGenerator.vech(rand(size(m_observed,1))),1);
            cvx_begin quiet
            variable v(s_n)
            minimize( s_alpha*vec((m_observed)*(m_observed)')'*m_duplication*v+s_beta*v'*(m_duplication')*m_duplication*v)
            subject to
            m_A*v==[k;zeros(size(m_observed,1),1)];
            m_B*v<=0;
            cvx_end
            m_laplacian=GraphLearningSmoothSignalGraphGenerator.my_ivech(v,m_duplication);
        end
        
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
        
        function v=getDiagIndices(n)
            v=zeros(n,1); 
            for i=1:n 
                if i==1 
                    v(i)=1; 
                else 
                    v(i)=i+v(i-1); 
                end 
            end 
        end
        
        function m_A=getLinEqualities(m_laplacian)
            X=GraphLearningSmoothSignalGraphGenerator.vech(m_laplacian);
            % %A contais the following info:
            % % tr(L)=n
            % % diag elements are in the positions prev_diag+diag_index
            v=GraphLearningSmoothSignalGraphGenerator.getDiagIndices(size(m_laplacian,1));
            indic=zeros(size(X,1),1);
            indic(v)=1;
            indic=indic';
            % %indic contains the tr(L)
            % % L*1=0
            
            m_A=zeros(size(m_laplacian,1)+1,size(X,1));
            v=v';
            for i=1:size(m_laplacian,1)
                if(i~=1)
                    ind=[v(i)-(0:i-1),v(i:end-1)+i];
                else
                    ind=horzcat(v(i),v(i:end-1)+i);
                end
                m_A(i+1,ind)=1;
            end
            m_A(1,:)=indic;
            %A*X
        end
        
        function m_B=getLinInequalities(m_laplacian)
            % B contains info about the Lij <= of zero
            % so must contain a line for each of these elements of X
            X=GraphLearningSmoothSignalGraphGenerator.vech(m_laplacian);
            v=GraphLearningSmoothSignalGraphGenerator.getDiagIndices(size(m_laplacian,1));
            m_B=zeros(size(X,1)-size(m_laplacian,1),size(X,1));
            for i=1:size(X,1)
                if ismember(i,v)
                    i=i-1;
                else
                    m_B(i,i)=1;
                end
            end
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
        
        function m_estimated=signalUpd(m_X,m_laplacian,s_alpha,m_W)
            % min norm(W*(X-Y),'fro')^2+a*tr(Y'*L*Y)
            % wrt Y
            % convex has closed form solution
            % m_W indicate positions of missing values
            m_estimated = zeros(size(m_X));
            if m_W == ones(size(m_W))
                m_estimated = (eye(size(m_laplacian))+s_alpha*m_laplacian)\m_X;
            else
                for m=1: size(m_estimated,2)
                    m_estimated(:,m) = (diag(m_W(:,1)) + ...
                        s_alpha*m_laplacian) \ (diag(m_W(:,1))*m_X(:,m));
                end
            end
        end
        
        function L=my_ivech(x,M)
            K2=size(x,1);
            K=(-1+sqrt(1+8*K2))/2;
            d=GraphLearningSmoothSignalGraphGenerator.getDiagIndices(K);
            vec=M*x;
            L=reshape(vec,K,K);
        end
            
        function m_dublication=dup(m_laplacian)
            % converts the vectorized form of lower triangular
            % part of L to the vectorized form of the full matrix
            % M*vech(L)=vec(L)
            
            X=GraphLearningSmoothSignalGraphGenerator.vech(m_laplacian);
            m_dublication=zeros(size(m_laplacian,1)^2,size(X,1));
            for i=1:size(m_laplacian,1)
                if i==1
                    v(i)=1;
                else
                    v(i)=i+v(i-1);
                end
            end
            for i=1:size(m_laplacian,1)
                if(i~=1)
                    ind=[v(i)-(0:i-1),v(i:end-1)+i];
                else
                    ind=horzcat(v(i),v(i:end-1)+i);
                end
                ind=sort(ind);
                for j=1:size(m_laplacian,1)
                    m_dublication(j+(i-1)*size(m_laplacian,1),ind(j))=1;
                end
            end
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
            % to the initial matrix
            %if size(v_stackedData,2)>size(v_stackedData,1)
            %    v_stackedData=v_stackedData';
            %end
            %if size(v_stackedData,2)~=1
            %    error('STACKED_DATA must be a column vector.')
            %end
            v_stackedData = v_stackedData(:);
            
            len = length(v_stackedData);
            N = (-1+sqrt(1+8*len))/2;
            
            if floor(N)~=N
                error(['The number of elemeents in STACKED_DATA must be conformable to' ...
                    'the inverse vech operation.'])
            end
            % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % %Input Checking
            % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % % %Initialize the output data
            m_matrixData=zeros(N);
            
            % % %Use a logical trick to inverse the vech
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
        
        function v_vec = vec(m_X)
            % contains the vectorized form of m_X
            % [m,n] = size(m_X);
            % v_vec = reshape(m_X,m*n,1);
            v_vec = m_X(:);
        end
    end
end
