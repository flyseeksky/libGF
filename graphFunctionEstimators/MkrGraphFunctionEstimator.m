classdef MkrGraphFunctionEstimator < GraphFunctionEstimator
    % Function estimator using multi-kernel regression method
    
    properties
        c_parsToPrint  = {'ch_name', 'legendString',};
		c_stringToPrint  = {'', ''};
		c_patternToPrint = {'%s%s', '%s%s'};
    end
    
    properties
		ch_name = 'Multi-kernel RR';
        m_kernel   %  N x N x P tensor, where 
		           %       N: number of vertices
				   %       P: number of kernels
        s_mu
        s_sigma    % only valid for single kernel
		s_type = 'RKHS superposition'; 
		           % can be 'RKHS superposition' or 'kernel superposition'
		
				   
    end
    
    properties (Dependent)
        legendString
    end
    
    methods
        function str = get.legendString(obj)
            % return string that can be used for legend generating
            % because when replicate along m_kernel, this property cannot 
            % be printed, thus define this dependent property to accomplish
            % this task
            nKernels = size(obj.m_kernel,3);
            if nKernels == 1
                str = sprintf('1 kernel, \\sigma = %3.2f', obj.s_sigma);
            else
                str = sprintf('%d kernels', nKernels);
            end
        end
    end
    
    methods
		
        function obj = MkrGraphFunctionEstimator(varargin)  % constructor
            obj@GraphFunctionEstimator(varargin{:});
        end
        
        function [m_estimate, m_alpha] = estimate(obj, m_samples, m_positions)
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
            [m_estimate, m_alpha] = obj.estimateSignal(m_samples, m_positions);
        end
        
        function [m_estimate, m_alpha] = estimateSignal(obj,m_samples,m_positions)
            % This is the function that does the real job
            [N,M,nKernel] = size(obj.m_kernel);  % N is # of vertices
            assert(N==M, 'Kernel matrix should be square');
			
            K = obj.m_kernel;
            K_observed = obj.m_kernel(m_positions, m_positions, :);
			S = size(m_positions, 1);
			
			kernelScale = NaN(nKernel,1);     % normalize kernel matrix
			for iKernel = 1 : size(K_observed,3)
				kernelScale(iKernel) = trace(K_observed(:,:,iKernel))/N;
				K_observed(:,:,iKernel) = K_observed(:,:,iKernel)/kernelScale(iKernel);
				K(:,:,iKernel) = K(:,:,iKernel) / kernelScale(iKernel);
			end
            
            % estimate alpha
            alpha = real( obj.estimateAlpha( m_samples, K_observed ) );
			
			% recover signal on whole graph
            m_estimate = zeros(N,1);  % y = \sum K_i * alpha_i
			m_alpha = zeros(length(m_positions), nKernel);
            for iKernel = 1 : nKernel
                alpha_i = alpha( ((iKernel-1)*S + 1) : iKernel*S ); % extract ai
                Ki = K(:,m_positions,iKernel); % kernel of the whole graph
                m_estimate = m_estimate + Ki*alpha_i;
                
				m_alpha(:,iKernel) = alpha_i;  % record alpha_i
            end
        end
        
        function a = estimateAlpha(obj, m_samples, K)
            % use group lasso to solve alpha
            % Input:
            %       m_samples       observed signal
            %       K               kernel matrix for observed part
            % Output:
            %       a               alpha in group lasso format
            %
            
            S = length(m_samples);
            u = obj.s_mu;             % regularization parameter
            y = m_samples;            % observed signal
            nKernel = size(K,3);      % # of kernels

            % change variable to group lasso format
            A = NaN(S,S*nKernel);
            for iKernel = 1 : nKernel
                Ki = K(:,:,iKernel);
                A(:, ((iKernel-1)*S + 1) : iKernel*S ) = mpower(Ki,1/2);
            end
% 			A = reshape(K, [size(K,1) size(K,2)*size(K,3)]);
            
            % set the parameter for group lasso solver
            lambda = S/2 * u;
            p = ones(nKernel,1) * S;
            rho = 1;
            alpha = 1;
            
            % solve the problem
            [r, history] = group_lasso(A, y, lambda, p, rho, alpha);
            r = real(r);
            
            % interpret the result
            a = NaN(S*nKernel,1);
            for iKernel = 1 : nKernel				
                sqrtKi = A(:, ((iKernel-1)*S+1) : iKernel*S);
                ri = r( ((iKernel-1)*S+1) : iKernel*S );
                a( ((iKernel-1)*S+1) : iKernel*S ) = sqrtKi \ ri;
            end
        end
        
	end
end