classdef RidgeRegressionGraphFunctionEstimator < GraphFunctionEstimator
    % Function estimator using kernel regression
    
    properties
        c_parsToPrint  = {'ch_name'};
		c_stringToPrint  = {''};
		c_patternToPrint = {'%s%s'};
    end
    
    properties
		ch_name = 'Kernel RR';
        m_kernel   %  N x N kernel matrix
		           %       N: number of vertices
				   
    end
        
    methods
		
        function obj = RidgeRegressionGraphFunctionEstimator(varargin)  % constructor
            obj@GraphFunctionEstimator(varargin{:});
        end
        
        function [v_estimate, m_alpha ] = estimate(obj, v_samples, v_positions)
			%    v_positions     S x 1 vector with the indices of the
			%                    sampled vertices
			% 
			%    v_samples:      S x 1 vector with the (possibly noisy)
			%                    function values: v_samples(i) is the
			%                    observation at node v_positions(i)
			%
			%
			%    v_estimate:     N x 1 vector with the signal estimate			
			%
			%    m_alpha:        S x 1 vector with the vector of alphas
			
			
			% Initial checks
            if isempty(obj.m_kernel) || isempty(obj.s_regularizationParameter)
                error('MkrGraphFunctionEstimator:notEnoughInfo',...
                    'Kernel and mu not set');
            elseif ~isequaln(size(v_samples),size(v_positions))				
                error('MkrGraphFunctionEstimator:inconsistentParameter',...
                    'size of m_positions and m_samples not the same');
            elseif max(v_positions(:)) > size(obj.m_kernel,1)
                error('MkrGraphFunctionEstimator:outOfBound', ...
                    'position out of bound');
            end						
            [N,Np] = size(obj.m_kernel);  % N is # of vertices
            assert(N==Np, 'Kernel matrix should be square');
			assert(size(v_samples,2)==1,'not implemented');
			assert(size(obj.m_kernel,3)==1);
			
			% Kernel preparation 	
            if length(obj.s_regularizationParameter)>1
				% choose s_regularizationParameter via cross validation
				s_mu = obj.crossValidation(v_samples,v_positions,obj.s_regularizationParameter);
			elseif length(obj.s_regularizationParameter) == 1
				s_mu = obj.s_regularizationParameter;
			else
				error('empty obj.s_regularizationParameter');
            end
            obj.s_regularizationParameter = s_mu;
			m_alpha = (obj.m_kernel(v_positions, v_positions) + size(v_samples,1)*s_mu*eye(size(v_samples,1)))\v_samples;
			v_estimate = obj.m_kernel(:,v_positions)*m_alpha;
			
			
			
		end
        
		function N = getNumOfVertices(obj)
            N = size(obj.m_kernel,1);
        end
		
	end
	
	
end