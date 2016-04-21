classdef MkrGraphFunctionEstimator < GraphFunctionEstimator
    % Function estimator using multi-kernel regression method
    
    properties
        c_parsToPrint  = {};
		c_stringToPrint  = {};
		c_patternToPrint = {};
    end
    
    properties
        m_kernel
        s_mu
    end
    
    methods
        function obj = MkrGraphFunctionEstimator(varargin)  % constructor
            obj@GraphFunctionEstimator(varargin{:});
        end
        
        function m_estimate = estimate(obj, m_samples, m_positions)
            if isempty(obj.m_kernel) || isempty(obj.s_mu)
                error('MkrGraphFunctionEstimator:notEnoughInfo',...
                    'Kernel and mu not set');
            elseif ~isequaln(size(m_samples),size(m_positions))
                error('MkrGraphFunctionEstimator:inconsistentParameter',...
                    'size of m_positions and m_samples not the same');
            elseif max(m_positions(:)) > size(obj.m_kernel,1)
                error('MkrGraphFunctionEstimator:outOfBound', ...
                    'position out of bound');
            end
            m_estimate = obj.estimateSignal(m_samples, m_positions);
        end
        
        function m_estimate = estimateSignal(obj,m_samples,m_positions)
            [N,~,nKernel] = size(obj.m_kernel);
            K_observed = obj.m_kernel(m_positions, m_positions, :);
            K = obj.m_kernel;
            S = size(m_positions, 1);
            % estimate a
            y = m_samples;
            a = obj.estimateAlpha( y, K_observed );
            %alpha_mat = reshape(a,S,nKernel);
            % get the estimated signal on the whole graph
            m_estimate = zeros(N,1);
            for iKernel = 1 : nKernel
                % extend ai from size S to size N
                alpha = zeros(N,1);
                ai = a( (iKernel - 1)*S + 1 : iKernel*S );
                alpha(m_positions) = ai;

                % get estimated signal
                Ki = K(:,:,iKernel);
                m_estimate = m_estimate + Ki*alpha;
            end
        end
        
        function a = estimateAlpha(obj, m_samples, K)   % estimate alpha
            S = length(m_samples);
            u = obj.s_mu;
            y = m_samples;
            nKernel = size(K,3);      % # of kernels

            % change variable to group lasso format
            A = nan(S,S*nKernel);
            for iKernel = 1 : nKernel
                Ki = K(:,:,iKernel);
                A(:, (iKernel-1)*S + 1 : iKernel*S ) = mpower(Ki,1/2);
            end
            % set the parameter for group lasso solver
            lambda = S/2 * u;
            p = ones(nKernel,1) * S;
            rho = 1;
            alpha = 1;
            % solve the problem
            [r, history] = group_lasso(A, y, lambda, p, rho, alpha);
            r = real(r);
            % interpret the result
            a = nan(S*nKernel,1);
            for iKernel = 1 : nKernel
                sqrtKi = A(:, (iKernel-1)*S + 1 : iKernel*S);
                ri = r( (iKernel-1)*S + 1 : iKernel*S );
                a( (iKernel-1)*S + 1 : iKernel*S ) = sqrtKi \ ri;
            end
        end
        
    end
    
    methods(Static)
        
        function m_kernel = getKernelFromFunction(m_laplacian,c_functions)
           % 
           % m_kernel       N x N x P tensor where m_kernel(:,:,p) is the
           %                kernel matrix for the graph laplacian kernel
           %                corresponding to the function of eigenvalues
           %                c_functions(p).
           %                N = size(m_laplacian，1）
           %                P = length(c_functions)
           %
                                
            
        end
        
        function m_kernel = getDiffusionKernel(m_laplacian,v_sigma2)
           % 
           % m_kernel       N x N x P tensor where m_kernel(:,:,p) is the
           %                kernel matrix for the graph diffusion kernel
           %                with parameter v_sigma2(p).
           %                N = size(m_laplacian，1）
           %                P = length(v_sigma2)
           %
           
           for k=length(v_sigma2):-1:1
               c_functions(k) = @(lambda) exp(v_sigma2(k)/(2*lambda));
           end
           
           m_kernel =  MkrGraphFunctionEstimator.getKernelFromFuncion(m_laplacian,c_functions);
           
        end
        
        
        
        
        
    end
    
end