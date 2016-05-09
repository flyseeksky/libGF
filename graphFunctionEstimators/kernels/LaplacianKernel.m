classdef LaplacianKernel < Parameter
    % KernelGenerator
    %       Class to generator graph kernel matrices
    %   
    
    properties
        c_parsToPrint = {'ch_type','s_sigma'};
		c_stringToPrint = {'','\sigma'};
		c_patternToPrint = {'%s kernel','%s = %f'};
    end
    
    properties
        h_r_inv         % cell array of r^{-1}(\lambda) in graph kernel
        m_laplacian     % laplacian matrix
                       
    end
    
    methods
        function obj = LaplacianKernel(varargin)
            obj@Parameter(varargin{:});
        end
        
        function t_kernelMatrix = getKernelMatrix(obj)
            % t_kernelMatrix is an N x N x P tensor where
            %        N = size(obj.m_laplacian,1)
            %        P = length(obj.h_r_inv)
            %
            
            if ~iscell(obj.h_r_inv)
                error('h_r_inv must be a cell array');
            end
            
            if isempty(obj.m_laplacian)
                error('Property m_laplacian cannot be empty');
            end
            
            [V,D] = eig(obj.m_laplacian);
            d = diag(D);  d(1) = 0;                  % fix a bug since the first
                                            % eigenvalue of L is always 0
            N = size(obj.m_laplacian,1);
            P = length(obj.h_r_inv);
            t_kernelMatrix = NaN(N,N,P);
            for p = 1 : P
                r = obj.h_r_inv{p};
                Kp = V * diag(r(d)) * V.';
                t_kernelMatrix(:,:,p) = ( Kp + Kp' )/2;
            end
        end
    end
    
    % type specific methods
    methods(Static)
        % diffusion kernel: r(lambda) = exp(sigma^2*lambda/2)
        function h_r_inv = diffusionKernelFunctionHandle(sigmaArray)
                       
            for iSigma = length(sigmaArray):-1:1
                sigma = sigmaArray(iSigma);
                h_r_inv{iSigma} = @(lambda) exp( - sigma^2 * lambda / 2 );
            end
            
        end
        
         % regularized kernel: r(lambda) = 1 + sigma^2*lambda
        function h_r_inv = regularizedKernelFunctionHandle( sigmaArray)
            
            for iSigma = length(sigmaArray):-1:1
                sigma = sigmaArray(iSigma);
                h_r_inv{iSigma} = @(lambda) 1./(1 + sigma^2 * lambda);
            end
            
        end
        
        
    end
    methods
        
       
    end
    
end