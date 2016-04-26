classdef KernelGenerator < Parameter
    % KernelGenerator
    %       Class to generator graph kernel matrices
    %   
    
    properties
        c_parsToPrint = {'ch_type'};
		c_stringToPrint = {''};
		c_patternToPrint = {'%s kernel'};
    end
    
    properties
        h_r             % cell array of r(\lambda) in graph kernel
        m_laplacian     % laplacian matrix
        ch_type         % type of kernel, can be diffusion, regularized, 
                        % one-step random walk or inverse cosine
    end
    
    methods
        function obj = KernelGenerator(varargin)
            obj@Parameter(varargin{:});
        end
        
        function t_kernelMatrix = getKernelMatrix(obj)
            if ~iscell(obj.h_r)
                error('h_r must be a cell array');
            end
            
            if isempty(obj.m_laplacian)
                error('Property m_laplacian cannot be empty');
            end
            
            [V,D] = eig(obj.m_laplacian);
            N = size(obj.m_laplacian,1);
            P = length(obj.h_r);
            t_kernelMatrix = NaN(N,N,P);
            for p = 1 : P
                r = obj.h_r{p};
                t_kernelMatrix(:,:,p) = V * r(D) * V';
            end
        end
    end
    
    % type specific methods
    methods
        % diffusion kernel: r(lambda) = exp(sigma^2*lambda/2)
        function t_kernelMatrix = getDiffusionKernel(obj, sigmaArray)
            if ~strcmp(obj.ch_type, 'diffusion')
                error('Kernel type not consistent');
            end
            
            for iSigma = 1 : length(sigmaArray)
                sigma = sigmaArray(iSigma);
                obj.h_r{iSigma} = @(lambda) exp( - sigma^2 * lambda / 2 );
                t_kernelMatrix = obj.getKernelMatrix();
            end
        end
        
        % regularized kernel: r(lambda) = 1 + sigma^2*lambda
        function t_kernelMatrix = getRegularizedKernel(obj, sigmaArray)
            if ~strcmp(obj.ch_type, 'regularized')
                error('Kernel type %s is not consistent with regularized',...
                    obj.ch_type);
            end
            
            for iSigma = 1 : length(sigmaArray)
                sigma = sigmaArray(iSigma);
                obj.h_r{iSigma} = @(lambda) 1./(1 + sigma^2 * lambda);
                t_kernelMatrix = obj.getKernelMatrix();
            end
        end
    end
    
end