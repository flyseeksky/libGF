%
%  FIGURES FOR THE PAPER ON MULTIKERNEL
%
% 

classdef MultikernelSimulations < simFunctionSet
	
	properties
	
	end
	
	methods
		
		% This is a very simple simulation to test bandlimited LS
		% estimation
		function F = compute_fig_1001(obj,niter)
			
			N = 100; % number of vertices
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.3,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',30);
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_numberOfSamples',40,'s_SNR',20);
			
			% 3. define graph function estimator
			estimator = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors',functionGenerator.basis);
			
			% Simulation
			m_graphFunction = functionGenerator.realization();
			[m_samples,m_positions] = sampler.sample(m_graphFunction);
			m_graphFunctionEstimate = estimator.estimate(m_samples,m_positions);
			
			% Performance assessment
			error = norm(m_graphFunctionEstimate - m_graphFunction,'fro')^2/size(m_graphFunction,1)
			
			F = F_figure('X',1:N,'Y',[m_graphFunctionEstimate,m_graphFunction]','leg',{'estimate','true'},'xlab','VERTEX','ylab','FUNCTION');
			
		end
		
		% This is a simple simulation to construct a Monte Carlo figure
		function F = compute_fig_2001(obj,niter)
			
			N = 100; % number of vertices
			S_vec = 10:10:100; % number of samples
			B = 20; % bandwidth of the estimated function
			B_vec = 10:10:50; % assumed bandwidth for estimation
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.3,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			bandlimitedFunctionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			graphFunction = bandlimitedFunctionGenerator.realization();
			generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',graphFunction);
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',20);
			sampler = sampler.replicate([],{},'s_numberOfSamples',num2cell(S_vec));
						
			% 3. define graph function estimator
			estimator = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors',bandlimitedFunctionGenerator.basis(N));
			estimator = estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});

			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,estimator);
			mse = Simulator.computeMse(res,Results('stat',graphFunction));			
			
			% Representation of results
			F = F_figure('X',Parameter.getXAxis(generator,sampler,estimator),'Y',mse,'leg',Parameter.getLegend(generator,sampler,estimator),'xlab',Parameter.getXLabel(generator,sampler,estimator),'ylab','MSE');
			
		end
		
		% This is a simple simulation to construct a Monte Carlo figure
		% Different from 2001, objets of different classes are concatenated
		function F = compute_fig_2002(obj,niter)
						
			N = 100; % number of vertices			
			B = 20; % bandwidth of the estimated function
			B_vec =         [10 20 30 10 20 30]; % assumed bandwidth for estimation
			SNR_vec = [15 25 15 15 15 25 25 25]; % SNR for each curve (first 2 for multikernel)
			
			S_vec = 10:10:100;
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.9,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			bandlimitedFunctionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			graphFunction = bandlimitedFunctionGenerator.realization();
			generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',graphFunction);			
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',20);			
			sampler = sampler.replicate('s_SNR',num2cell(SNR_vec),'s_numberOfSamples',num2cell(S_vec));		
						
			% 3. BL graph function estimator
			bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors',bandlimitedFunctionGenerator.basis(N));			
			bl_estimator.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator = bl_estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});
					
			% 4. MKL function estimator
		    m_laplacian = bandlimitedFunctionGenerator.basis(N);
			m_kernel = cat(3,pinv(m_laplacian)+1e-10*eye(N),pinv(m_laplacian^2)+1e-10*eye(N));
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',m_kernel,'s_mu',1e-5);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};

			est = [mkl_estimator;mkl_estimator;bl_estimator];
			
			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,est);
			mse = Simulator.computeMse(res,Results('stat',graphFunction));
			
			% Representation			
			F = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est));
			
		end
		
		% Figure to illustrate the interpolating functions (columns of the
		% kernel matrix)
		function F = compute_fig_2003(obj,niter)
			
			vertexNum = 100;
			columnInd = 50;
			
			sigma2 = .1;
			rDiffusionKernel = @(lambda,sigma2) exp(sigma2*lambda/2);
			KcolDiffusionKernel = MultikernelSimulations.columnLaplacianKernelCircularGraph(vertexNum,@(lambda) rDiffusionKernel(lambda,sigma2) , columnInd);
			
			% computation via analytic expression
			epsilon = 1e-6;
			rLaplacianReg = @(lambda,epsilon) lambda + epsilon;
			KcolLaplacianReg_analytic = MultikernelSimulations.columnLaplacianKernelCircularGraph(vertexNum,@(lambda) rLaplacianReg(lambda,epsilon) , columnInd);
			
			% direct computation
			A = circshift(eye(vertexNum),1)+circshift(eye(vertexNum),-1);
			L = diag(sum(A,2))-A;
			h_rFun = @(lambda) rLaplacianReg(lambda,epsilon);
			kG = KernelGenerator('m_laplacian',L,'h_r',{h_rFun});
			m_KernelMatrix = inv(kG.getKernelMatrix);
			KcolLaplacianReg_direct = m_KernelMatrix(:,columnInd);
			
						
			%F = F_figure('X',1:vertexNum,'Y',[KcolLaplacianReg_analytic';KcolLaplacianReg_direct']);			
			multiplot_array(1) = F_figure('X',1:vertexNum,'Y',[KcolLaplacianReg_direct']);
			multiplot_array(2) = F_figure('X',1:vertexNum,'Y',[KcolLaplacianReg_analytic']);
			F = F_figure('multiplot_array',multiplot_array);
		end
		
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  simulations with MKL on synthetic data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% Simulation to test the regularization parameter
		function F = compute_fig_3001(obj,niter)
            % 
			
			SNR = 20; % dB
			N = 100;
            S_Vec = 10:20:80;
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.2,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',30);
			m_graphFunction = functionGenerator.realization();
            generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
			% 3. generate Kernel matrix
			sigmaArray = linspace(0.01, 1.5, 30);
			L = graph.getLaplacian();
            kG = KernelGenerator('ch_type','diffusion','m_laplacian',L);
			m_kernel = kG.getDiffusionKernel(sigmaArray);
			
			% 4. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);
            sampler = sampler.replicate('s_numberOfSamples', num2cell(S_Vec),[],{}); 
			
			% 5. define function estimator
            estimator = MkrGraphFunctionEstimator('s_mu',1e-3);
            estimator = estimator.replicate([],{}, ...
                'm_kernel', mat2cell(m_kernel, N, N, ones(1,size(m_kernel,3))));
			
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter, m_graphFunction);
            
            % Representation
            F = F_figure('X',sigmaArray,'Y',mse, ...
                'leg',Parameter.getLegend(generator,sampler, estimator),...
                'xlab','\sigma','ylab','Normalized MSE');		  
        end	
        
        function F = compute_fig_3002(obj, niter)
            SNR = 20; % dB
			N = 100;
            u_Vec = logspace(-6,0,50);
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.1,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',30);
			m_graphFunction = functionGenerator.realization();
            generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
			% 3. generate Kernel matrix
			sigmaArray = linspace(0.01, 1.5, 20);
            %sigmaArray = 0.80;
			L = graph.getLaplacian();
            kG = KernelGenerator('ch_type','diffusion','m_laplacian',L);
			m_kernel = kG.getDiffusionKernel(sigmaArray);
            
            % 4. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples',50);
            
            % 5. define function estimator
            estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel);
            estimator = estimator.replicate([],{}, ...
                's_mu', num2cell(u_Vec));
			
            [m_samples, m_positions] = sampler.sample(m_graphFunction);
			m_alpha = zeros( length(m_samples), size(m_kernel,3), length(u_Vec) );
			for i = 1 : length(u_Vec)
				estimator_now = estimator(i);
				
				[~, alpha] = estimator_now.estimate(m_samples, m_positions);
				m_alpha(:,:,i) = alpha;
			end
			
			anorm = sum( m_alpha.^2, 1 );
			anorm = permute(anorm, [3 2 1]);
            
            for i = 1:length(sigmaArray)
                legendStr{i} = sprintf('\\sigma=%2.2f',sigmaArray(i));
            end
			
			F = F_figure('X', u_Vec, 'Y', anorm', 'logx', true, ...
				'xlab', '\mu', 'ylab', '||\alpha_i||^2','leg',legendStr);

		end
		
		function F = compute_fig_3003(obj, niter)
			% find best regularization
            SNR = 20; % dB
			N = 100;
            u_Vec = logspace(-6,0,50);
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.1,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',30);
			m_graphFunction = functionGenerator.realization();
            generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
			% 3. generate Kernel matrix
			sigmaArray = linspace(0.1, 1.5, 20);
            %sigmaArray = 0.80;
			L = graph.getLaplacian();
            kG = KernelGenerator('ch_type','diffusion','m_laplacian',L);
			m_kernel = kG.getDiffusionKernel(sigmaArray);
            
            % 4. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples',40);
            
            % 5. define function estimator
            estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel);
            estimator = estimator.replicate([],{}, ...
                's_mu', num2cell(u_Vec));
			
			
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter, m_graphFunction);
			
			F = F_figure('X', u_Vec, 'Y', mse, 'logx', true, ...
				'xlab', '\mu', 'ylab', 'MSE');
        end
        
        function F = compute_fig_4001(obj,niter)
            
			
			SNR = 20; % dB
			N = 100;
            S_Vec =  10:10:80;
            v_bandwidth = [2 5 10 20 40];
            mu_Vec = [1e-2 1e-2 1e-2 0.003];
            
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.1,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',30);
			m_graphFunction = functionGenerator.realization();
            generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
            L = graph.getLaplacian();
            
			% 4. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);
            sampler = sampler.replicate([],{}, 's_numberOfSamples', num2cell(S_Vec)); 
			
			% 5. define function estimator
            bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors', L);
            bl_estimator = bl_estimator.replicate('s_bandwidth', ...
                num2cell(v_bandwidth), [], {});
            
            % 3. generate Kernel matrix
			
            kG = KernelGenerator('ch_type','diffusion','m_laplacian',L);
            sigmaArray = [0.86 0.80 0 0];
            c_kernel{1} = kG.getDiffusionKernel(sigmaArray(1));
            c_kernel{2} = kG.getDiffusionKernel(sigmaArray(2));
            sigmaArray2 = [3 0.8];
            c_kernel{3} = kG.getDiffusionKernel(sigmaArray2);
            sigmaArray20 = linspace(0.1,1.5,20); %[0.1 0.3 0.5 0.8 0.95 1.1 1.3 1.5];
			c_kernel{4} = kG.getDiffusionKernel(sigmaArray20);
            
            for i = 1 : length(sigmaArray)
                mk_estimator(i) = MkrGraphFunctionEstimator('s_mu',mu_Vec(i),...
                    's_sigma',sigmaArray(i), 'm_kernel', c_kernel{i}, ...
                    'c_replicatedVerticallyAlong', {'legendString'});
            end
            
            %estimator = [bl_estimator; mk_estimator(:)];
            estimator = mk_estimator(:);
		
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter, m_graphFunction);
            
            % Representation
            F = F_figure('X',S_Vec,'Y',mse, ...
                'leg',Parameter.getLegend(generator,sampler, estimator),...
                'xlab','sample size','ylab','Normalized MSE');	  
        end
		
		
		
		% old sim. to compare with 4003
		function F = compute_fig_4002(obj,niter)
			
			% regularization parameter
			u = 1e-6;
			N = 100;
% 			p = 0.4;
% 			atilde = (1:N)';
% 			alpha = exp( - atilde / 50 );
% 			[V,d] = generateGraphData(N,p);
						
% 			N = size(V,1);
			bandwidthArray = 30:10:60;
			nBandwidth = length(bandwidthArray);
			sampleSizeArray = 40; %20:2:90;
			nSampleSize = length(sampleSizeArray);
			MONTE_CARLO = niter;
			
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.8,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',30);
			m_graphFunction = functionGenerator.realization();
			ftrue = m_graphFunction;
			L = graph.getLaplacianEigenvectors();
			[V,D] = eig(L);
			d = diag(D);
			
			% generate signal
			% alpha = zeros(N,1) ; %1./(1:N)'; %sort(1+2*rand(B,1),'descend');
			% alpha(1:4) = 1;
			% atilde = (1:N)';
			% alpha = exp( - atilde / 50 );
			% [V,d] = generateGraphData(N,p);
			% ftrue = 10*V*alpha;
			SNR = 20; % dB
			
			% generate Kernel matrix
			sigma1Array = linspace(0.01, 0.5, 10);
			sigma2Array = linspace(0.01, 1.5, 3);
			nSigma1 = length(sigma1Array);
			nSigma2 = length(sigma2Array);
			K1 = zeros(N,N,nSigma1);
			K2 = zeros(N,N,nSigma2);
			for iSigma = 1 : nSigma1
				sigma = sigma1Array( iSigma );
				K1(:,:,iSigma) = V * diag( exp( - sigma^2/2 * d ) ) * V';
			end
			for iSigma = 1 : nSigma2
				sigma = sigma2Array( iSigma );
				K2(:,:,iSigma) = V * diag( exp( - sigma^2/2 * d ) ) * V';
			end
			
			
			
			%NMSE_MK = nan(nSampleSize);
% 			NMSE_LS = nan(nSampleSize, nBandwidth);
			for iSampleSize = 1:nSampleSize
				S = sampleSizeArray (iSampleSize);
				%
				% simulate MKL estimator
				%
				estimator1 = @(sampleSet, observedSignal) MKL_estimator(sampleSet, ...
					observedSignal, K1, u);
				NMSE_MK1(iSampleSize)  = ...
					MultikernelSimulations.sim_MKL(ftrue,S, SNR, estimator1, MONTE_CARLO);
				estimator2 = @(sampleSet, observedSignal) MKL_estimator(sampleSet, ...
					observedSignal, K2, u);
				NMSE_MK2(iSampleSize)  = ...
					MultikernelSimulations.sim_MKL(ftrue,S, SNR, estimator2, MONTE_CARLO);
				
				for iBandwidth = 1 : nBandwidth
					B = bandwidthArray( iBandwidth );
					
					fprintf('Progress: %3.1f%%\n', ...
						100*( (iSampleSize-1)*nBandwidth + iBandwidth ) / ...
						( nBandwidth*nSampleSize) );
					%
					% least square estimator
					%
% 					if S < B
% 						continue;
% 					end
% 					estimator = @(sampleSet, signal) LS_estimator(V,B,sampleSet,signal);
% 					NMSE_LS(iSampleSize,iBandwidth) = ...
% 						MultikernelSimulations.sim_MKL(ftrue,S,SNR,estimator,MONTE_CARLO);
				end
				
			end
			
			NMSE_MK1
			NMSE_MK2
% 			plot(sampleSizeArray, NMSE_LS, 'linewidth',1);
% 			hold on
% 			plot(sampleSizeArray, NMSE_MK1, 'linewidth',2);
% 			hold on
% 			plot(sampleSizeArray, NMSE_MK2, 'linewidth',2);
% 			hold off
% 			xlabel('sample size')
% 			ylabel('NMSE')
% 			ylim([0 1])
% 			for iBandwidth = 1:nBandwidth
% 				legendString{iBandwidth} = sprintf('B=%d',bandwidthArray(iBandwidth));
% 			end
% 			legendString{end+1} = 'MK method with 10 kernels';
% 			legendString{end+1} = 'MK method with 3 kernels';
% 			legend(legendString);
			
			
			F = [];
		end
		
		
		% new sim
        function F = compute_fig_4003(obj,niter)
            
			
			SNR = 20; % dB
			N = 100;
            S_Vec = 40; % 10:10:80;
            v_bandwidth = [2 5 10 20 40];
            mu = 5e-6;
            
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.8,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',30);
			m_graphFunction = functionGenerator.realization();
            generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
            L = graph.getLaplacian();
            
			% 4. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);
            sampler = sampler.replicate([],{}, 's_numberOfSamples', num2cell(S_Vec)); 
			
			% 5. define function estimator
            bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors', L);
            bl_estimator = bl_estimator.replicate('s_bandwidth', ...
                num2cell(v_bandwidth), [], {});
            
            % 3. generate Kernel matrix
			
            kG = KernelGenerator('ch_type','diffusion','m_laplacian',L);
            sigmaArray = [0.01 0.80 0 0];
            c_kernel{1} = kG.getDiffusionKernel(sigmaArray(1));
            c_kernel{2} = kG.getDiffusionKernel(sigmaArray(2));
            sigmaArray2 = linspace(0.01,1.5, 2);
            c_kernel{3} = kG.getDiffusionKernel(sigmaArray2);
            sigmaArray20 = linspace(0.01, 1.5, 20);
			c_kernel{4} = kG.getDiffusionKernel(sigmaArray20);
            
            for i = 1:4
                mk_estimator(i,1) = MkrGraphFunctionEstimator('s_mu',mu,...
                    's_sigma',sigmaArray(i), 'm_kernel', c_kernel{i}, ...
                    'c_replicatedVerticallyAlong', {'legendString'});
            end
            
            %estimator = [bl_estimator; mk_estimator(:)];
            estimator = mk_estimator;
		
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter, m_graphFunction)
            
            % Representation
%             F = F_figure('X',S_Vec,'Y',mse, ...
%                 'leg',Parameter.getLegend(generator,sampler, estimator),...
%                 'xlab','sample size','ylab','Normalized MSE');	  
           F = [];
        end
		
		
		
		
	end
	
	
	methods(Static)
		
		% =========================================================================
		% utility functions
		% =========================================================================
		function NMSE = sim_MKL(trueSignal,S, SNR,estimator,MONTE_CARLO)
			signalPower = norm(trueSignal)^2/length(trueSignal);
			noisePower = signalPower / 10^(SNR/10);
			
			N = length(trueSignal);
			N_SE = zeros(MONTE_CARLO,1);
			for iMonteCarlo = 1 : MONTE_CARLO
				% random generate a sample set
				componentArray = partition_set(N, S);
				sampleSet = componentArray(1).index;
				
				% generate observed signal
				observedSignal = trueSignal(sampleSet) + ...
					sqrt(noisePower) * randn(S,1);
				
				% estimate signal using the estimator
				estimatedSignal = estimator( sampleSet, observedSignal );
				
				% compute square error
				N_SE(iMonteCarlo) = norm(estimatedSignal - trueSignal)^2 / norm(trueSignal)^2;
			end
			
			NMSE = median(N_SE);
			
		end
		
		function Kcol = columnLaplacianKernelCircularGraph(vertexNum,rFun,columnInd)
			% Kcol is a vertexNum x 1 vector that corresponds to the
			% columnInd-th column of the Laplacian kernel matrix of a
			% circular graph when the r function is rFun. 
			%
			% rFun must accept vector-valued inputs.
			
			Dinds = (1:vertexNum)-columnInd;
			
			for rowInd = 1:vertexNum				
				Kcol(rowInd,1) = (1/vertexNum)*sum( exp(1j*2*pi/vertexNum*(0:vertexNum-1)*Dinds(rowInd))./rFun(2*(1-cos(2*pi/vertexNum*(0:vertexNum-1)))));
			end
			Kcol = real(Kcol);
		end
		
	end

	
end