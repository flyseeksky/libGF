%
%  FIGURES FOR THE PAPER ON MULTIKERNEL
%
%  TSP paper figures: 1006, 3101, 3305
%

classdef MultikernelSimulations < simFunctionSet
	
	properties
	
	end
	
	methods
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  1. Generic simulations
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% This is a very simple simulation to test bandlimited LS
		% estimation
		function F = compute_fig_1001(obj,niter)
			
			N = 100; % number of vertices
			B = 30;  % bandwidth
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.3,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_numberOfSamples',40,'s_SNR',20);
			
			% 3. define graph function estimator
			estimator = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',B);
			
			% Simulation
			m_graphFunction = functionGenerator.realization();
			[m_samples,m_positions] = sampler.sample(m_graphFunction);
			m_graphFunctionEstimate = estimator.estimate(m_samples,m_positions);
			
			% Performance assessment
			error = norm(m_graphFunctionEstimate - m_graphFunction,'fro')^2/size(m_graphFunction,1)
			
			F = F_figure('X',1:N,'Y',[m_graphFunctionEstimate,m_graphFunction]','leg',{'estimate','true'},'xlab','VERTEX','ylab','FUNCTION');
			
		end
					
		% This is a very simple simulation to test the computation of the
		% cut-off frequency in [narang2013structured] and [anis2016proxies]
		function F = compute_fig_1002(obj,niter)
			
			N = 100; % number of vertices
			B = 30;  % bandwidth
			SNR = 10; % dB
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.3,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_numberOfSamples',40,'s_SNR',SNR);
			
			% 3. define graph function estimator
			estimator_known_freq = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',B);
			estimator_unknown_freq = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',-1);
			
			% Simulation
			m_graphFunction = functionGenerator.realization();
			[m_samples,m_positions] = sampler.sample(m_graphFunction);
			m_graphFunctionEstimate_known_freq = estimator_known_freq.estimate(m_samples,m_positions);
			m_graphFunctionEstimate_unknown_freq = estimator_unknown_freq.estimate(m_samples,m_positions);
			
			% Performance assessment
			error_known_freq = norm(m_graphFunctionEstimate_known_freq - m_graphFunction,'fro')^2/size(m_graphFunction,1)
			error_unknown_freq = norm(m_graphFunctionEstimate_unknown_freq - m_graphFunction,'fro')^2/size(m_graphFunction,1)
			
			F = F_figure('X',1:N,'Y',[m_graphFunction,m_graphFunctionEstimate_known_freq,m_graphFunctionEstimate_unknown_freq]','leg',{'true','estimate (known freq.)','estimate (unknown freq.)'},'xlab','VERTEX','ylab','FUNCTION','styles',{'-','--','-.'});
			
		end
				
		% This is a simple simulation to construct a Monte Carlo figure
		function F = compute_fig_1003(obj,niter)
			
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
			estimator = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);
			estimator = estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});

			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,estimator);
			mse = Simulator.computeNmse(res,Results('stat',graphFunction));			
			
			% Representation of results
			F = F_figure('X',Parameter.getXAxis(generator,sampler,estimator),'Y',mse,'leg',Parameter.getLegend(generator,sampler,estimator),'xlab',Parameter.getXLabel(generator,sampler,estimator),'ylab','MSE','ylimit',[0 1.5]);
			
		end
		
		% This is a simple simulation to construct a Monte Carlo figure
		% Different from 2001, objets of different classes are concatenated
		function F = compute_fig_1004(obj,niter)
						
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
			bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);			
			bl_estimator.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator = bl_estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});
					
			% 4. MKL function estimator
		    m_laplacian = bandlimitedFunctionGenerator.basis(N);
			m_kernel = cat(3,pinv(m_laplacian)+1e-10*eye(N),pinv(m_laplacian^2)+1e-10*eye(N));
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',m_kernel,'s_regularizationParameter',1e-5);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};

			est = [mkl_estimator;mkl_estimator;bl_estimator];
			
			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,est);
			mse = Simulator.computeNmse(res,Results('stat',graphFunction));
			
			% Representation			
			F = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est));
			
		end
		
		% Figure to check analytic expression for interpolating functions
		% (columns of the kernel matrix) in a circular graph
		function F = compute_fig_1005(obj,niter)
			
			vertexNum = 100;
			columnInd = 25;
			A = circshift(eye(vertexNum),1)+circshift(eye(vertexNum),-1);
			L = diag(sum(A,2))-A;
						
			% Computation through analytic expression for
			% a) Laplacian regularization
			epsilon = .01;
			rLaplacianReg = @(lambda,epsilon) lambda + epsilon;
			KcolLaplacianReg_analytic = MultikernelSimulations.columnLaplacianKernelCircularGraph(vertexNum,@(lambda) rLaplacianReg(lambda,epsilon) , columnInd);
			
			% b) Diffusion kernel
			sigma2 = 3;
			rDiffusionKernel = @(lambda,sigma2) exp(sigma2*lambda/2);
			KcolDiffusionKernel_analytic = MultikernelSimulations.columnLaplacianKernelCircularGraph(vertexNum,@(lambda) rDiffusionKernel(lambda,sigma2) , columnInd);
						
			% direct computation for
			
			% a) regularized laplacian
			h_rFun_inv = @(lambda) 1./rLaplacianReg(lambda,epsilon);
			kG = LaplacianKernel('m_laplacian',L,'h_r_inv',{h_rFun_inv});
			m_KernelMatrix = kG.getKernelMatrix;
			KcolLaplacianReg_direct = m_KernelMatrix(:,columnInd);
			
			% a) diffusion kernel
			h_rFun_inv = @(lambda) 1./rDiffusionKernel(lambda,sigma2);
			kG = LaplacianKernel('m_laplacian',L,'h_r_inv',{h_rFun_inv});
			m_KernelMatrix = kG.getKernelMatrix;
			KcolDiffusionKernel_direct = m_KernelMatrix(:,columnInd);
			
			
			F(1) = F_figure('X',1:vertexNum,'Y',[KcolLaplacianReg_direct';KcolLaplacianReg_analytic'],'styles',{'-','--'});
			F(2) = F_figure('X',1:vertexNum,'Y',[KcolDiffusionKernel_direct';KcolDiffusionKernel_analytic'],'styles',{'-','--'});
			
		end
		
		% Figure to illustrate the interpolating functions (columns of the
		% kernel matrix) in a circular graph
		function F = compute_fig_1006(obj,niter)
			
			vertexNum = 100;
			columnInd = 25;
			A = circshift(eye(vertexNum),1)+circshift(eye(vertexNum),-1);
			L = diag(sum(A,2))-A;
						
			% Computation through analytic expression for
			% a) Laplacian regularization
			rLaplacianReg = @(lambda,s2) 1+s2*lambda;
			v_sigma2_LaplacianReg = [1 20 100];
			for i_sigma2 = length(v_sigma2_LaplacianReg):-1:1				
				KcolLaplacianReg(i_sigma2,:) = MultikernelSimulations.columnLaplacianKernelCircularGraph(vertexNum,@(lambda) rLaplacianReg(lambda,v_sigma2_LaplacianReg(i_sigma2)) , columnInd)';
				leg{i_sigma2} = sprintf('Laplacian reg. (\\sigma^2 = %g )',v_sigma2_LaplacianReg(i_sigma2));
			end
			KcolLaplacianReg = diag(1./max(KcolLaplacianReg,[],2))*KcolLaplacianReg;
			
			
			% b) Diffusion kernel			
			v_sigma2_DiffusionKernel = [1 20 100];
			rDiffusionKernel = @(lambda,sigma2) exp(sigma2*lambda/2);
			i_legLen = length(leg);
			for i_sigma2 = length(v_sigma2_DiffusionKernel):-1:1
				KcolDiffusionKernel(i_sigma2,:) = MultikernelSimulations.columnLaplacianKernelCircularGraph(vertexNum,@(lambda) rDiffusionKernel(lambda,v_sigma2_DiffusionKernel(i_sigma2)) , columnInd)';
				leg{i_sigma2+i_legLen} = sprintf('Diffusion k. (\\sigma^2 = %g )',v_sigma2_DiffusionKernel(i_sigma2));
			end		
			KcolDiffusionKernel = diag(1./max(KcolDiffusionKernel,[],2))*KcolDiffusionKernel;
			
			caption = sprintf('%d-th column of the kernel matrix for a circular graph with N = %d vertices.',columnInd,vertexNum);
			m_Y = [KcolLaplacianReg;KcolDiffusionKernel];
			F = F_figure('X',1:2:vertexNum,'Y',m_Y(:,1:2:vertexNum),'leg',leg,'styles',{'-','-x','-o','--','--x','--o'},'colorp',3,'xlab','Vertex index (n)','ylab','Function value','caption',caption,'leg_pos_vec',[0.5546    0.5271    0.2333    0.3715]);
			
		end
		
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  2. simulations with estimators for bandlimited signals on
		% %%  synthetic data 
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
		% Simple simulation to test [narang2013structured] and
		% [anis2016proxies] cut-off freq. estimation method
		function F = compute_fig_2001(obj,niter)
						
			N = 100; % number of vertices			
			B = 20; % bandwidth of the estimated function
			B_vec =         [20]; % assumed bandwidth for estimation
			SNR_vec = [15 25 25 25]; % SNR for each curve (first 2 for multikernel)
			
			S_vec = 10:10:100;
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.9,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			m_laplacian = graph.getLaplacian(); 
			bandlimitedFunctionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			graphFunction = bandlimitedFunctionGenerator.realization();
			generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',graphFunction);			
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',20);			
			sampler = sampler.replicate('s_SNR',num2cell(SNR_vec),'s_numberOfSamples',num2cell(S_vec));		
						
			% 3. BL graph function estimator
			bl_estimator_known_freq = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);			
			bl_estimator_known_freq.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator_known_freq = bl_estimator_known_freq.replicate('s_bandwidth',num2cell(B_vec),'',{});
					
			% 4. BL estimator with unknown frequency
			bl_estimator_unknown_freq = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',-1);			
			bl_estimator_unknown_freq.c_replicatedVerticallyAlong = {'ch_name','s_bandwidth'};
						
			% 5. MKL function estimator		    
			sigma2Array = linspace(0.1, .5 , 20);            
            kG = LaplacianKernel('m_laplacian',m_laplacian,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigma2Array));
			m_kernel = kG.getKernelMatrix();
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',m_kernel,'s_regularizationParameter',1e-3);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};

			est = [mkl_estimator;mkl_estimator;bl_estimator_known_freq;bl_estimator_unknown_freq];
			
			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,est);
			mse = Simulator.computeNmse(res,Results('stat',graphFunction));
			
			% Representation			
			F = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est),'ylimit',...
				[0 1.5],'ylab','NMSE','tit',Parameter.getTitle(graphGenerator,bandlimitedFunctionGenerator,generator,sampler));
			
		end
				
		
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% %%  3. simulations with MKL on synthetic data
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
        % Index
        %       sigma w.r.t bandwidth              Figure 3100
        %       sigma w.r.t sample size            Figure 3101
        %       performance comparison of
        %       single kernel and multikernel      Figure 3102
        %       IIA with bandlimited kernels       Figure 3103
        %       RKHS with bandlimited kernels      Figure 3104
        %       Sparsity of alpha w.r.t mu         Figure 3201
        %       NMSE vs mu(regularization)         Figure 3202
        %       Test parameter for Cortes          Figure 3203
        
        
        % 1) Figures for tuning kernel parameters==========================
        
        % Figure: NMSE vs sigma
        %    Instead of drawing one curve per sample size, draw one curve
        %    per bandwidth (of signal).
        function F = compute_fig_3100(obj,niter)		
			[N,p,SNR,sampleSize,~] = MultikernelSimulations.simulationSetting();
            mu = 1e-4;
            bandwidthVec = [5 10 20 30 40];
						
			% generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            
            % generate signal on this graph
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph);
            functionGenerator.b_generateSameFunction = 1;
            generator = functionGenerator.replicate('s_bandwidth', ...
                num2cell(bandwidthVec), [], {} );
			
			% generate Kernel matrix
			sigmaArray = sqrt(linspace(0.01, 1.5, 30));
            %sigma = 0.8;
			L = graph.getLaplacian();
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray));
			m_kernel = kG.getKernelMatrix();
            
			% define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, ...
                's_numberOfSamples', sampleSize);
			
			% define function estimator
            estimator = MkrGraphFunctionEstimator('s_regularizationParameter',mu);
            estimator = estimator.replicate([],{}, ...
                'm_kernel', mat2cell(m_kernel, N, N, ones(1,size(m_kernel,3))));
			
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter);
            
            % Representation
            F = F_figure('X',sigmaArray.^2,'Y',mse, ...
                'leg',Parameter.getLegend(generator,sampler, estimator),...
                'xlab','\sigma^2','ylab','Normalized MSE',...
                'tit', sprintf('N=%d, p=%2.2f, \\mu=%3.1d, S = %d', N, p, mu, sampleSize));		  
		end
		
		% Figure NMSE vs samplesize
		%    each curve corresponds to one bandwidth -- bandlimited kernel
		function F = compute_fig_3110(obj,niter)
			[N,p,SNR,sampleSize,~] = MultikernelSimulations.simulationSetting();
			mu = 1e-4;
			B = 20;  % used to generate graph signals
            bandwidthVec = 10:5:40;  % used to build bandlimited kernels
			S_vec = 10:5:90;
			beta = 1000;
						
			% generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            
            % generate signal on this graph
			generator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth', B);
			
			% generate Kernel matrix
			L = graph.getLaplacian();
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.bandlimitedKernelFunctionHandle( L , bandwidthVec , beta));
			m_kernel = kG.getKernelMatrix();
            
			% define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);
			sampler = sampler.replicate([],{},'s_numberOfSamples', num2cell(S_vec));
			
			% define function estimator
            estimator = MkrGraphFunctionEstimator('s_regularizationParameter',mu);
            estimator = estimator.replicate('m_kernel', mat2cell(m_kernel, N, N, ones(1,size(m_kernel,3))),[],{});
			
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter);
            
            % Representation
			for i = 1:length(bandwidthVec)
				legendStr{i}  = sprintf('B = %d', bandwidthVec(i));
			end
            F = F_figure('X',S_vec,'Y',mse, ...
                'leg',legendStr, 'ylimit', [0 1.5], ...
                'xlab','sample size','ylab','Normalized MSE',...
                'tit', sprintf('N=%d, p=%2.2f, \\mu=%3.1d, signal bandwidth = %d', N, p, mu, B));		  
		end
		
		% Figure: ||alpha_i|| vs \mu
		%      shows sparsity pattern of bandlimited kernels
		%      This simulation aims to show that the best single kernel is the
		%      to vanish when mu is increased large enough. If this is true, 
		%	   then bandwidth estimation is possible
		function F = compute_fig_3130(obj, niter)
			
			[N,p,SNR,sampleSize,~] = MultikernelSimulations.simulationSetting();
			B = 20;  % used to generate graph signals
			u_Vec = logspace(-6,0,50);
			
			s_beta = 1e3;				% amplitude of the bandlimited kernel
			v_B_values = 10:5:40;		% bandwidth for the bandlimited kernel
			
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
			
			% 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			v_graphFunction = functionGenerator.realization();
			
			% 3. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples',sampleSize);
			
			% 4. generate Kernel matrix
			kG = LaplacianKernel('m_laplacian',graph.getLaplacian(),'h_r_inv',LaplacianKernel.bandlimitedKernelFunctionHandle(graph.getLaplacian(),v_B_values,s_beta));
			m_kernel = kG.getKernelMatrix();
			
			% 5. define function estimator
			estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel);
			estimator = estimator.replicate([],{}, 's_regularizationParameter', num2cell(u_Vec));
			
			[m_samples, m_positions] = sampler.sample(v_graphFunction);
			m_alpha = zeros( length(m_samples), size(m_kernel,3), length(u_Vec) );
			for icount = 1 : length(u_Vec)
				estimator_now = estimator(icount);
				
				[v_graphFunction_estimate, alpha] = estimator_now.estimate(m_samples, m_positions);
				m_alpha(:,:,icount) = alpha;
				
				v_nmse(icount) = norm( v_graphFunction - v_graphFunction_estimate)^2/norm( v_graphFunction )^2;
			end
			
			anorm = sum( m_alpha.^2, 1 );
			anorm = permute(anorm, [3 2 1]);			
			 
            for icount = 1:length(v_B_values)
                legendStr{icount} = sprintf('B = %d',v_B_values(icount));
            end
			
			multiplot_array(1,1) = F_figure('X', u_Vec, 'Y', anorm', 'logx', true, ...
				'xlab', '\mu', 'ylab', '||\alpha_i||^2','leg',legendStr,'leg_pos','East',...
				'tit', sprintf('signal B = %d, S = %d',B, sampleSize));
			multiplot_array(2,1) = F_figure('X', u_Vec, 'Y', v_nmse, 'logx', true, ...
				'xlab', '\mu', 'ylab', 'NMSE');
			F(1) = F_figure('multiplot_array',multiplot_array);
						
			F(2) = F_figure('Y',graph.getFourierTransform(v_graphFunction)','tit','Fourier transform of target signal','xlab','Freq. index','ylab','Function value');
			
			
		end
				
		% Figure: NMSE vs sigma (diffusion kernel parameter)
		% This figure will show the importance of choosing the right
		%   parameter (sigma for diffusion kernel, may change to other
		%     parameter if different kernel types are used.)
		function F = compute_fig_3101(obj,niter)
			[N,p,SNR,sampleSize,bandwidth] = MultikernelSimulations.simulationSetting();
            S_Vec = 10:10:80;
            mu = 1e-2;
						
			% generate graph and signal
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
			%functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',30);
			functionGenerator = ExponentiallyDecayingGraphFunctionGenerator('graph',graph,'s_bandwidth',30,'s_decayingRate',.5);
			m_graphFunction = functionGenerator.realization();
            generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
			% 3. generate Kernel matrix
			sigmaArray = sqrt(linspace(0.01, 1, 30));
			L = graph.getLaplacian();
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray));
			m_kernel = kG.getKernelMatrix();
            
            
			% 4. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);
            sampler = sampler.replicate('s_numberOfSamples', num2cell(S_Vec),[],{}); 
			
			% 5. define function estimator
            estimator = MkrGraphFunctionEstimator('s_regularizationParameter',mu);
            estimator = estimator.replicate([],{}, ...
                'm_kernel', mat2cell(m_kernel, N, N, ones(1,size(m_kernel,3))));
			
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter);
            
            % Representation
            F = F_figure('X',sigmaArray.^2,'Y',mse, ...
                'leg',Parameter.getLegend(generator,sampler, estimator),...
                'xlab','\sigma^2','ylab','Normalized MSE',...
                'tit', sprintf('N=%d, p=%2.2f, \\mu=%3.1d', N, p, mu),...
				'leg_pos','northwest');		  
		end			
		
		% Figure: NMSE vs S (single kernel with different paramter and multi-kernel
        %         with different number of kernels
        function F = compute_fig_3102(obj,niter)
            [N,p,SNR,~,bandwidth] = MultikernelSimulations.simulationSetting();                 
                                             % from Figure 3100 the best kernel
                                             % for bandwidht = 30 is sigma = 0.8
            mu_Vec = 1e-2*ones(5,1);
            S_Vec = 10:10:80;
            
						
			% generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            
            % generate graph function
			generator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',bandwidth);
% 			m_graphFunction = functionGenerator.realization();
%             generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
 			
            L = graph.getLaplacian();
            
			% define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);
            sampler = sampler.replicate([],{}, 's_numberOfSamples', num2cell(S_Vec)); 
			

			% kernels for single kernel estimators
            sigmaArray = sqrt([0.2 0.80 1 0 0]);
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray(1)));			
            c_kernel{1} = kG.getKernelMatrix();
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray(2)));			
            c_kernel{2} = kG.getKernelMatrix();    
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray(3)));			
            c_kernel{3} = kG.getKernelMatrix();
            
			% kernels for multi-kernel estimators
            sigmaArray2 = sqrt([0.4 0.8 1.2]);
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray2));			
            c_kernel{4} = kG.getKernelMatrix();
            
            sigmaArray20 = sqrt(linspace(0.1,1.5,20)); %[0.1 0.3 0.5 0.8 0.95 1.1 1.3 1.5];
			kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray20));			
            c_kernel{5} = kG.getKernelMatrix();
            
            %c_kernel{4} = kG.getDiffusionKernel(sigmaArray20);
            
            for i = 1 : length(sigmaArray)
                estimator(i,1) = MkrGraphFunctionEstimator('s_regularizationParameter',mu_Vec(i),...
                    's_sigma',sigmaArray(i), 'm_kernel', c_kernel{i}, ...
                    'c_replicatedVerticallyAlong', {'legendString'});
            end
            
            
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter);
            tit = Parameter.getTitle(graphGenerator,generator,sampler, estimator);
			
            % Representation
            F = F_figure('X',S_Vec,'Y',mse, ...
                'leg',Parameter.getLegend(generator,sampler, estimator),...
                'xlab','sample size','ylab','Normalized MSE',...
                'tit', tit);	  
		end
		
		% Simulation to see how IIA works with bandlimited kernels
		% Figure: |theta_i| vs mu for i = 1,..,#kernels
		% Depicts the pattern  of theta in IIA
		% as the regularization paramter mu increases
		function F = compute_fig_3103(obj, niter)
			
            SNR = 20; % dB
			N = 100;
			B = 30; % bandwidth
			S = 80; % number of observed vertices
            u_Vec = logspace(-6,6,50);
			
			s_beta = 10000; % amplitude parameter of the bandlimited kernel
			v_B_values = 10:10:50; % bandwidth parameter for the bandlimited kernel
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.5,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
			
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			v_graphFunction = functionGenerator.realization();
            %generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
			% 3. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples',S);
			
			% 4. generate Kernel matrix
            kG = LaplacianKernel('m_laplacian',graph.getLaplacian(),'h_r_inv',LaplacianKernel.bandlimitedKernelFunctionHandle(graph.getLaplacian(),v_B_values,s_beta));
			m_kernel = kG.getKernelMatrix();                   
            
            % 5. define function estimator
            estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel,'ch_type','kernel superposition');
            estimator = estimator.replicate([],{}, ...
                's_regularizationParameter', num2cell(u_Vec));
			
            [m_samples, m_positions] = sampler.sample(v_graphFunction);
			m_theta = zeros( size(m_kernel,3), length(u_Vec) );
			v_nmse = zeros( 1 , length(u_Vec) );
			for icount = 1 : length(u_Vec)				
				[v_graphFunction_now,~,m_theta(:,icount)] = estimator(icount).estimate(m_samples, m_positions);				 
				v_nmse(icount) = norm( v_graphFunction - v_graphFunction_now)^2/norm( v_graphFunction )^2;
			end
			
            
            for icount = 1:length(v_B_values)
                legendStr{icount} = sprintf('B = %2.2f',v_B_values(icount));
            end
			
			multiplot_array(1,1) = F_figure('X', u_Vec, 'Y', m_theta, 'logx', true, ...
				'xlab', '\mu', 'ylab', 'Entries of \theta','leg',legendStr,'leg_pos','West');
			multiplot_array(2,1) = F_figure('X', u_Vec, 'Y', v_nmse, 'logx', true, ...
				'xlab', '\mu', 'ylab', 'NMSE');
			F(1) = F_figure('multiplot_array',multiplot_array);
			
			F(2) = F_figure('Y',graph.getFourierTransform(v_graphFunction)','tit','Fourier transform of target signal','xlab','Freq. index','ylab','Function value');
			

		end

		% Simulation to see how MKL with 'RKHS superposition' works with
		% bandlimited kernels  
		% Figure:  ||alpha_i|| vs mu for i = 1,..,#kernels
		% Depicts the sparsity pattern  of theta in IIA
		% as regularization paramter mu increases, theta would become more
		% more sparse
		function F = compute_fig_3104(obj, niter)
			
			SNR = 20; % dB
			N = 100;
			B = 20; % bandwidth
			S = 50; % number of observed vertices
			u_Vec = logspace(-6,0,50);
			
			s_beta = 1e6; % amplitude parameter of the bandlimited kernel
			v_B_values = 10:5:30; % bandwidth parameter for the bandlimited kernel
			
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.5,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
			
			% 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			v_graphFunction = functionGenerator.realization();
			%generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
			% 3. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples',S);
			
			% 4. generate Kernel matrix
			kG = LaplacianKernel('m_laplacian',graph.getLaplacian(),'h_r_inv',LaplacianKernel.bandlimitedKernelFunctionHandle(graph.getLaplacian(),v_B_values,s_beta));
			m_kernel = kG.getKernelMatrix();
			
			% 5. define function estimator
			estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel);
			estimator = estimator.replicate([],{}, ...
				's_regularizationParameter', num2cell(u_Vec));
			
			[m_samples, m_positions] = sampler.sample(v_graphFunction);
			m_alpha = zeros( length(m_samples), size(m_kernel,3), length(u_Vec) );
			for icount = 1 : length(u_Vec)
				estimator_now = estimator(icount);
				
				[v_graphFunction_estimate, alpha] = estimator_now.estimate(m_samples, m_positions);
				m_alpha(:,:,icount) = alpha;
				
				v_nmse(icount) = norm( v_graphFunction - v_graphFunction_estimate)^2/norm( v_graphFunction )^2;
			end
			
			anorm = sum( m_alpha.^2, 1 );
			anorm = permute(anorm, [3 2 1]);			
			 
            for icount = 1:length(v_B_values)
                legendStr{icount} = sprintf('B = %d',v_B_values(icount));
            end
			
			multiplot_array(1,1) = F_figure('X', u_Vec, 'Y', anorm', 'logx', true, ...
				'xlab', '\mu', 'ylab', '||\alpha_i||^2','leg',legendStr,'leg_pos','East');
			multiplot_array(2,1) = F_figure('X', u_Vec, 'Y', v_nmse, 'logx', true, ...
				'xlab', '\mu', 'ylab', 'NMSE');
			F(1) = F_figure('multiplot_array',multiplot_array);
						
			F(2) = F_figure('Y',graph.getFourierTransform(v_graphFunction)','tit','Fourier transform of target signal','xlab','Freq. index','ylab','Function value');
			
			
		end

		
		% print version of 3100
		function F = compute_fig_3105(obj, niter)
			F = obj.load_F_structure(3100);
			F.styles = {'-','--','-^','--^','-*'};
			F.xlab = '\sigma^2';
			F.caption = F.tit;
			F.ylab = 'NMSE';
			F.tit = '';
			F.leg_pos_vec = [ 0.3073    0.6041    0.2206    0.3045];
        end
        
        % 
        function F = compute_fig_3120(obj,niter)		
			[N,p,SNR,sampleSize,~] = MultikernelSimulations.simulationSetting();
            % mu = 1e-4;
            sampleSize = 40;
            bandwidthVec = [5 10 20 30 40];
						
			% generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            
            % generate signal on this graph
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph);
            functionGenerator.b_generateSameFunction = 1;
            generator = functionGenerator.replicate('s_bandwidth', ...
                num2cell(bandwidthVec), [], {} );
			
			% generate Kernel matrix
			%sigmaArray = sqrt(linspace(0.01, 1.5, 10));
			sigmaArray = sqrt(linspace(0.1, 1.5, 20));
            %sigma = 0.8;
			L = graph.getLaplacian();
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray));
			m_kernel = kG.getKernelMatrix();
            
			% define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples', sampleSize);
            	
			% define function estimator
            estimatorTemp = MkrGraphFunctionEstimator('s_numFoldValidation',10);
            estimator = estimatorTemp.replicate([],{}, ...
                'm_kernel', mat2cell(m_kernel, N, N, ones(1,size(m_kernel,3))));
            
			% cross validation
            v_mu = logspace(-10,0,14);
			m_graphFunction = generator(1).realization();
            for estIndex = 1 : length(estimator)
                [v_samples, v_positions] = sampler(1).sample(m_graphFunction);
                mu(estIndex) = estimator(estIndex).crossValidation(v_samples, v_positions, v_mu);
                estimator(estIndex).s_regularizationParameter = mu(estIndex);
            end
            
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter);
            
            % Representation
            F = F_figure('X',sigmaArray.^2,'Y',mse, ...
                'leg',Parameter.getLegend(generator,sampler, estimator),...
                'xlab','\sigma^2','ylab','Normalized MSE',...
                'tit', sprintf('N=%d, p=%2.2f, S = %d', N, p, sampleSize));		  
		end
		
		% 2) Figures for tuning the regularization parameter ==============
		
		% Figure: ||alpha_i|| vs mu
		% Depicts the sparsity pattern  of alpha
		% as regularization paramter mu increases, alpha would become more
		% more sparse, so more and more ||alpha_i|| will go to zero
		function F = compute_fig_3201(obj, niter)
			
            [N,p,SNR,sampleSize,bandwidth] = MultikernelSimulations.simulationSetting();
            u_Vec = logspace(-6,0,50);
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',bandwidth);
			m_graphFunction = functionGenerator.realization();
            generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
		
			% 3. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples',sampleSize);
			
			% 4. generate Kernel matrix
			sigmaArray = sqrt(linspace(0.01, 1.5, 20));
            %sigmaArray = 0.80;
			L = graph.getLaplacian();
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray));
			m_kernel = kG.getKernelMatrix();
            
            % 5. define function estimator
            estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel);
            estimator = estimator.replicate([],{}, ...
                's_regularizationParameter', num2cell(u_Vec));
			
            [m_samples, m_positions] = sampler.sample(m_graphFunction);
			m_alpha = zeros( length(m_samples), size(m_kernel,3), length(u_Vec) );
			for i = 1 : length(u_Vec)
				estimator_now = estimator(i);
				
				[~, alpha] = estimator_now.estimate(m_samples, m_positions);
				m_alpha(:,:,i) = alpha;
			end
			
			alphanorm = sum( m_alpha.^2, 1 );
			alphanorm = permute(alphanorm, [3 2 1]);
            
            for i = 1:length(sigmaArray)
                legendStr{i} = sprintf('\\sigma=%2.2f',sigmaArray(i));
            end
			
			F = F_figure('X', u_Vec, 'Y', alphanorm', 'logx', true, ...
				'xlab', '\mu', 'ylab', '||\alpha_i||^2','leg',legendStr,'leg_pos','West',...
                'tit',sprintf('N=%d,p=%2.2f,B=%d,S=%d',N,p,bandwidth,sampleSize));

		end
						
		% Figure: NMSE vs mu (regularization parameter)
		% Find the best regularization paramter for each method
		%    To find the best regularization paramter for other methods,
		%    only need to replace the estimator
		function F = compute_fig_3202(obj, niter)
						
            [N,p,SNR,sampleSize,bandwidth] = MultikernelSimulations.simulationSetting();
            u_Vec = logspace(-6,0,50);
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',bandwidth);
			m_graphFunction = functionGenerator.realization();
            generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
			% 3. generate Kernel matrix
			sigmaArray = linspace(0.1, 1.5, 20);
            %sigmaArray = 0.80;
			L = graph.getLaplacian();
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray));
			m_kernel = kG.getKernelMatrix();
            
            % 4. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples',sampleSize);
            
            % 5. define function estimator
            estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel);
            estimator = estimator.replicate([],{}, ...
                's_regularizationParameter', num2cell(u_Vec));
			
			
			% Simulation
            mse = Simulate(generator, sampler, estimator, niter);
			
			F = F_figure('X', u_Vec, 'Y', mse, 'logx', true, ...
				'xlab', '\mu', 'ylab', 'MSE', ...
                'tit', sprintf('N=%d, p=%2.2f, B=%d,S=%d, numOfKernels=%d', ...
                N, p, bandwidth,sampleSize, length(sigmaArray)));
		end
        
		% Simulation to test parameters for Cortes' MKL
		% Figure: |theta_i| vs mu for i = 1,..,#kernels
		% Depicts the pattern  of theta in IIA
		% as regularization paramter mu increases
		function F = compute_fig_3203(obj, niter)
			
            SNR = 20; % dB
			N = 100;
			B = 30; % bandwidth
            u_Vec = logspace(-6,6,50);
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.5,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
			
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			v_graphFunction = functionGenerator.realization();
            %generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
			% 3. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples',50);
			
			% 4. generate Kernel matrix
			sigmaArray = linspace(0.01, 1.5, 20);            			
            kG = LaplacianKernel('m_laplacian',graph.getLaplacian(),'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigmaArray));
			m_kernel = kG.getKernelMatrix();                   
            
            % 5. define function estimator
            estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel,'ch_type','kernel superposition');
            estimator = estimator.replicate([],{}, ...
                's_regularizationParameter', num2cell(u_Vec));
			
            [m_samples, m_positions] = sampler.sample(v_graphFunction);
			m_theta = zeros( size(m_kernel,3), length(u_Vec) );
			v_nmse = zeros( 1 , length(u_Vec) );
			for icount = 1 : length(u_Vec)				
				[v_graphFunction_now,~,m_theta(:,icount)] = estimator(icount).estimate(m_samples, m_positions);				 
				v_nmse(icount) = norm( v_graphFunction - v_graphFunction_now)^2/norm( v_graphFunction )^2;
			end
			
            
            for icount = 1:length(sigmaArray)
                legendStr{icount} = sprintf('\\sigma=%2.2f',sigmaArray(icount));
            end
			
			multiplot_array(1,1) = F_figure('X', u_Vec, 'Y', m_theta, 'logx', true, ...
				'xlab', '\mu', 'ylab', 'Entries of \theta','leg',legendStr,'leg_pos','West');
			multiplot_array(2,1) = F_figure('X', u_Vec, 'Y', v_nmse, 'logx', true, ...
				'xlab', '\mu', 'ylab', 'NMSE');
			F = F_figure('multiplot_array',multiplot_array);

        end
        
        %
        % Estimate bandwidth of bandlimited graph signals
        %     Use mkl estimator to find the best kernel, which is generated using
        %     the corresponding bandlimited kernels
        function F = compute_fig_3232(obj, niter)
						
            [N,p,SNR,sampleSize,bandwidth] = MultikernelSimulations.simulationSetting();
            u_Vec = logspace(-6,0,50);
            
						
			% 1. generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
            % 2. generate graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',bandwidth);
			%m_graphFunction = functionGenerator.realization();
            %generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',m_graphFunction);
			
			% 3. generate Kernel matrix
			%sigmaArray = linspace(0.1, 1.5, 20);
            B_vec = 10:10:90;
            beta = 10;   % for bandlimited kernel
            %sigmaArray = 0.80;
			L = graph.getLaplacian();
            kG = LaplacianKernel('m_laplacian',L,'h_r_inv',LaplacianKernel.bandlimitedKernelFunctionHandle(L, B_vec, beta));
			m_kernel = kG.getKernelMatrix();
            
            % 4. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR, 's_numberOfSamples',sampleSize);
            
            % 5. define function estimator
            estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel, 's_regularizationParameter', 1e-1);
            estimator.b_estimateFreq = 1;
            %estimator = estimator.replicate([],{}, ...
            %    's_regularizationParameter', num2cell(u_Vec));
			
			
			% Simulation
            bandwidth_vec = 10:10:60;
            estimated_bandwidth = zeros(length(bandwidth_vec),niter);
            for iBand = 1:length(bandwidth_vec)
                for iter = 1:niter
                    functionGenerator.s_bandwidth = bandwidth_vec(iBand);
                    v_graphFunction = functionGenerator.realization();
                    [m_samples, m_positions] = sampler.sample(v_graphFunction);
                    [v_estimate, m_alpha , v_theta, main_kernel_ind] = estimator.estimate(m_samples, m_positions);
                    estimated_bandwidth(iBand, iter) = B_vec(main_kernel_ind);
                end
                MultikernelSimulations.printSimulationProgress(iBand, iter, length(bandwidth_vec), niter)
            end
            est_bandwidth = mean(estimated_bandwidth, 2);
            %mse = Simulate(generator, sampler, estimator, niter);
			
			F = F_figure('X', bandwidth_vec, 'Y', [bandwidth_vec; est_bandwidth'], ...
				'xlab', 'experiment index', 'ylab', 'bandwidth', ...
                'tit', sprintf('N=%d, p=%2.2f, S=%d, numOfKernels=%d', ...
                N, p, sampleSize, length(B_vec)), ...
                'leg',{'true bandwidth','estimated bandwidth'});
		end
        
        

 	
		% 3) Figures to compare MKL and bandlimited =======================
		
		% Simple MC simulation to test MKL methods and compare them with
		% bandlimited estimators
		% - bandlimited signal, but MC does not average across signal
		%   realizations
		function F = compute_fig_3301(obj,niter)
						
			N = 100; % number of vertices			
			B = 20; % bandwidth of the estimated function
			B_vec =         [20]; % assumed bandwidth for estimation
			SNR_vec = [25 25 25 25]; % SNR for each curve (first 2 for multikernel)
			
			S_vec = 10:10:100;
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.9,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			m_laplacian = graph.getLaplacian(); 
			bandlimitedFunctionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			graphFunction = bandlimitedFunctionGenerator.realization();
			generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',graphFunction);			
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',20);			
			sampler = sampler.replicate('s_SNR',num2cell(SNR_vec),'s_numberOfSamples',num2cell(S_vec));		
						
			% 3. BL graph function estimator
			bl_estimator_known_freq = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);			
			bl_estimator_known_freq.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator_known_freq = bl_estimator_known_freq.replicate('s_bandwidth',num2cell(B_vec),'',{});
					
			% 4. BL estimator with unknown frequency
			bl_estimator_unknown_freq = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian,'s_bandwidth',-1);			
			bl_estimator_unknown_freq.c_replicatedVerticallyAlong = {'ch_name','s_bandwidth'};
						
			% 5. MKL function estimators		    
			sigma2Array = linspace(0.1, .5 , 20);            
            kG = LaplacianKernel('m_laplacian',m_laplacian,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigma2Array));
			m_kernel = kG.getKernelMatrix();
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',m_kernel,'s_regularizationParameter',1e-3);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};			
			mkl_estimator = mkl_estimator.replicate('ch_type',{'RKHS superposition','kernel superposition'},'',[]);

			est = [mkl_estimator;bl_estimator_known_freq;bl_estimator_unknown_freq];
			
			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,est);
			mse = Simulator.computeNmse(res,Results('stat',graphFunction));
			
			% Representation			
			F = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est),'ylimit',...
				[0 1.5],'ylab','NMSE','tit',Parameter.getTitle(graphGenerator,bandlimitedFunctionGenerator,generator,sampler));
			
		end
				
		% MC simulation to compare MKL and bandlimited estimators 
		% - bandlimited signal, but MC does not average across signal
		%   realizations
		function F = compute_fig_3302(obj,niter)
						
			N = 100; % number of vertices			
			B = 20; % bandwidth of the true function
			B_vec =         [10 20 30 -1]; % assumed bandwidth for estimation
			SNR = 5; % dB
			
			S_vec = 10:10:100;
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.7 ,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			m_laplacian = graph.getLaplacian(); 
			bandlimitedFunctionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			graphFunction = bandlimitedFunctionGenerator.realization();
			generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',graphFunction);			
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);			
			sampler = sampler.replicate('',[],'s_numberOfSamples',num2cell(S_vec));		
						
			% 3. BL graph function estimator
			bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);			
			bl_estimator.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator = bl_estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});
								
			% 4. MKL function estimators		    
			sigma2Array = linspace(0.1, .5 , 20);            
            kG = LaplacianKernel('m_laplacian',m_laplacian,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigma2Array));
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',kG.getKernelMatrix(),'s_regularizationParameter',5e-3);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};			
			mkl_estimator = mkl_estimator.replicate('ch_type',{'RKHS superposition','kernel superposition'},'',[]);

			est = [mkl_estimator;bl_estimator];
			
			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,est);
			mse = Simulator.computeNmse(res,Results('stat',graphFunction));
			
			% Representation			
			F = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est),'ylimit',...
				[0 1.5],'ylab','NMSE','tit',Parameter.getTitle(graphGenerator,bandlimitedFunctionGenerator,generator,sampler));
			
		end
		
		% MC simulation to compare MKL and bandlimited estimators
		% - bandlimited signal, but MC does AVERAGE across signal
		%   realizations (INCOMPLETE)
		function F = compute_fig_3303(obj,niter)
						
			N = 100; % number of vertices			
			B = 20; % bandwidth of the true function
			B_vec =         [10 20 30 -1]; % assumed bandwidth for estimation
			SNR = 5; % dB
			
			S_vec = 10:10:100;
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.7 ,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			m_laplacian = graph.getLaplacian(); 
			generator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);			
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);			
			sampler = sampler.replicate('',[],'s_numberOfSamples',num2cell(S_vec));		
						
			% 3. BL graph function estimator
			bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);			
			bl_estimator.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator = bl_estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});
								
			% 4. MKL function estimators		    
			sigma2Array = linspace(0.1, .5 , 20);            
            kG = LaplacianKernel('m_laplacian',m_laplacian,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigma2Array));
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',kG.getKernelMatrix(),'s_regularizationParameter',5e-3);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};			
			mkl_estimator = mkl_estimator.replicate('ch_type',{'RKHS superposition','kernel superposition'},'',[]);

			est = [mkl_estimator;bl_estimator];
			
			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,est);
			mse = Simulator.computeNmse(res,Results('stat',graphFunction));
			
			% Representation			
			F = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est),'ylimit',...
				[0 1.5],'ylab','NMSE','tit',Parameter.getTitle(graphGenerator,generator,sampler));
			
		end

		% MC simulation to compare MKL and bandlimited estimators 
		% - signal with exp. decaying spectrum. MC does not average across
		%   signal realizations
		function F = compute_fig_3304(obj,niter)
						
			N = 100; % number of vertices			
			B = 20; % bandwidth of the true function
			B_vec =         [10 20 30 -1]; % assumed bandwidth for estimation
			SNR = 5; % dB
			s_decayingRate = .5; % for decaying spectrum
			
			S_vec = 10:10:100;
			
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.7 ,'s_numberOfVertices',N);
			graph = graphGenerator.realization;
			m_laplacian = graph.getLaplacian(); 
			bandlimitedFunctionGenerator = ExponentiallyDecayingGraphFunctionGenerator('graph',graph,'s_bandwidth',B,'s_decayingRate',s_decayingRate);
			graphFunction = bandlimitedFunctionGenerator.realization();
			generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',graphFunction);			
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);			
			sampler = sampler.replicate('',[],'s_numberOfSamples',num2cell(S_vec));		
						
			% 3. BL graph function estimator
			bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);			
			bl_estimator.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator = bl_estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});
								
			% 4. MKL function estimators		    
			sigma2Array = linspace(0.1, .5 , 20);            
            kG = LaplacianKernel('m_laplacian',m_laplacian,'h_r_inv',LaplacianKernel.diffusionKernelFunctionHandle(sigma2Array));
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',kG.getKernelMatrix(),'s_regularizationParameter',5e-3);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};			
			mkl_estimator = mkl_estimator.replicate('ch_type',{'RKHS superposition','kernel superposition'},'',[]);

			est = [mkl_estimator;bl_estimator];
			
			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,est);
			mse = Simulator.computeNmse(res,Results('stat',graphFunction));
			
			% Representation			
			F(1) = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est),'ylimit',...
				[0 1.5],'ylab','NMSE','tit',Parameter.getTitle(graphGenerator,bandlimitedFunctionGenerator,generator,sampler),'leg_pos','southwest');
			F(2) = F_figure('Y',graph.getFourierTransform(graphFunction)','tit','Fourier transform of target signal','xlab','Freq. index','ylab','Function value');
			
		end
		
		
		% 4) Figures illustrate bandlimited kernels =======================
		
		% MC simulation to compare BL estimators and MKL estimators with BL
		% kernels. 
		% - bandlimited signal, but MC does not average across signal
		%   realizations
		function F = compute_fig_3401(obj,niter)
						
			N = 100; % number of vertices			
			B = 20; % bandwidth of the true function
			B_vec = [10 20 30 -1]; % assumed bandwidth for estimation
			SNR = 10; % dB
			S_vec = 10:5:100;
			
			s_beta = 1e5; % amplitude parameter of the bandlimited kernel
			v_B_values = 10:5:30; % bandwidth parameter for the bandlimited kernel
						
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.7 ,'s_numberOfVertices',N);
			graph = graphGenerator.realization;			
			bandlimitedFunctionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			graphFunction = bandlimitedFunctionGenerator.realization();
			generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',graphFunction);			
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);			
			sampler = sampler.replicate('',[],'s_numberOfSamples',num2cell(S_vec));		
						
			% 3. BL graph function estimator
			bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);			
			bl_estimator.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator = bl_estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});
								
			% 4. MKL function estimators	
			% 4. generate Kernel matrix
			kG = LaplacianKernel('m_laplacian',graph.getLaplacian(),'h_r_inv',LaplacianKernel.bandlimitedKernelFunctionHandle(graph.getLaplacian(),v_B_values,s_beta));			
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',kG.getKernelMatrix(),'s_regularizationParameter',5e-3);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};			
			mkl_estimator = mkl_estimator.replicate('ch_type',{'RKHS superposition','kernel superposition'},'',[]);
			est = [mkl_estimator;bl_estimator];
			
			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,est);
			mse = Simulator.computeNmse(res,Results('stat',graphFunction));
			
			% Representation			
			F(1) = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est),'ylimit',...
				[0 1.5],'ylab','NMSE','tit',Parameter.getTitle(graphGenerator,bandlimitedFunctionGenerator,generator,sampler));
			F(2) = F_figure('Y',graph.getFourierTransform(graphFunction)','tit','Fourier transform of target signal','xlab','Freq. index','ylab','Function value');
			
		end
		
		
		% MC simulation to compare BL estimators and MKL estimators with BL
		% kernels. 
		% - bandlimited signal, but MC does not average across signal
		%   realizations
		function F = compute_fig_3402(obj,niter)
						
			N = 100; % number of vertices			
			B = 20; % bandwidth of the true function
			B_vec = [10 20 30 -1]; % assumed bandwidth for estimation
			SNR = 10; % dB
			S_vec = 10:5:100;
			
			s_beta = 1e10; % amplitude parameter of the bandlimited kernel
			v_B_values = 10:5:30; % bandwidth parameter for the bandlimited kernel
						
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.7 ,'s_numberOfVertices',N);
			graph = graphGenerator.realization;			
			bandlimitedFunctionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			graphFunction = bandlimitedFunctionGenerator.realization();
			generator =  FixedGraphFunctionGenerator('graph',graph,'graphFunction',graphFunction);			
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);			
			sampler = sampler.replicate('',[],'s_numberOfSamples',num2cell(S_vec));		
						
			% 3. BL graph function estimator
			bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);			
			bl_estimator.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator = bl_estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});
								
			% 4. MKL function estimators	
			% 4. generate Kernel matrix
			kG = LaplacianKernel('m_laplacian',graph.getLaplacian(),'h_r_inv',LaplacianKernel.bandlimitedKernelFunctionHandle(graph.getLaplacian(),v_B_values,s_beta));			
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',kG.getKernelMatrix(),'s_regularizationParameter',5e-3);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};			
			mkl_estimator = mkl_estimator.replicate('ch_type',{'RKHS superposition','RKHS superposition','kernel superposition'},'',[]);
            mkl_estimator(1).b_finishSingleKernel = 1;
			mkl_estimator(1).s_finishRegularizationParameter = 1e7;
			mkl_estimator(1).c_replicatedVerticallyAlong = [mkl_estimator(1).c_replicatedVerticallyAlong,'b_finishSingleKernel'];
			est = [mkl_estimator;bl_estimator];
			
			% Simulation
			res = Simulator.simStatistic(niter,generator,sampler,est);
			mse = Simulator.computeNmse(res,Results('stat',graphFunction));
			
			% Representation			
			F(1) = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est),'ylimit',...
				[0 1.5],'ylab','NMSE','tit',Parameter.getTitle(graphGenerator,bandlimitedFunctionGenerator,generator,sampler));
			F(2) = F_figure('Y',graph.getFourierTransform(graphFunction)','tit','Fourier transform of target signal','xlab','Freq. index','ylab','Function value');
			
		end
		
		
		% MC simulation to compare BL estimators and MKL estimators with BL
		% kernels. 
		% - bandlimited signal, but MC does not average across signal
		%   realizations
		% - single kernel finish implements CV		
		function F = compute_fig_3403(obj,niter)
            [N,p,SNR,sampleSize, bandwidth] = MultikernelSimulations.simulationSetting();
            B = 20;
			B_vec = [10 20 30 -1]; % assumed bandwidth for estimation
			S_vec = 10:5:100;      
			
			s_beta = 1e5; % amplitude parameter of the bandlimited kernel
			v_B_values = 10:5:30; % bandwidth parameter for the bandlimited kernel
						
			% 1. define graph function generator
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', p ,'s_numberOfVertices',N);
			graph = graphGenerator.realization;			
			generator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',B);
			
			% 2. define graph function sampler
			sampler = UniformGraphFunctionSampler('s_SNR',SNR);			
			sampler = sampler.replicate('',[],'s_numberOfSamples',num2cell(S_vec));		
						
			% 3. BL graph function estimator
			bl_estimator = BandlimitedGraphFunctionEstimator('m_laplacian',graph.getLaplacian);			
			bl_estimator.c_replicatedVerticallyAlong = {'ch_name'};
			bl_estimator = bl_estimator.replicate('s_bandwidth',num2cell(B_vec),'',{});
								
			% 4. MKL function estimators	
			% 4. generate Kernel matrix
			kG = LaplacianKernel('m_laplacian',graph.getLaplacian(),'h_r_inv',LaplacianKernel.bandlimitedKernelFunctionHandle(graph.getLaplacian(),v_B_values,s_beta));			
			mkl_estimator = MkrGraphFunctionEstimator('m_kernel',kG.getKernelMatrix(),'s_regularizationParameter',1e-2);
			mkl_estimator.c_replicatedVerticallyAlong = {'ch_name'};			
			mkl_estimator = mkl_estimator.replicate('ch_type',{'RKHS superposition','RKHS superposition','kernel superposition'},'',[]);
            v_regPar = 10.^(-6:0);
			mkl_estimator(1).singleKernelPostEstimator= RidgeRegressionGraphFunctionEstimator('s_regularizationParameter',v_regPar);			
			mkl_estimator(1).c_replicatedVerticallyAlong = [mkl_estimator(1).c_replicatedVerticallyAlong,'s_regularizationParameter'];
			est = [mkl_estimator;bl_estimator];
			
			% Simulation
			%res = Simulator.simStatistic(niter,generator,sampler,est);
			%mse = Simulator.computeNmse(res,Results('stat',graphFunction));
            mse = Simulate(generator, sampler, est, niter);
			
			% Representation			
			F(1) = F_figure('X',Parameter.getXAxis(generator,sampler,est),...
                'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),...
                'xlab',Parameter.getXLabel(generator,sampler,est),'ylimit', [0 1.5],...
				'ylab','NMSE','tit',Parameter.getTitle(graphGenerator,generator,sampler));
			%F(2) = F_figure('Y',graph.getFourierTransform(graphFunction)','tit','Fourier transform of target signal','xlab','Freq. index','ylab','Function value');
			
		end
		
		
		
		
	end
	
	
	methods(Static)
		
		% =========================================================================
		% utility functions
		% =========================================================================
        function [N,p,SNR,sampleSize, bandwidth] = simulationSetting()
            % generate some commonly used paramters across simulations in order
            % to make all the simulations consistent
            N = 100;  % # of vertices
            p = 0.25; % edge existence prob of Erdos Renyi random graph model
            sampleSize = 40;
            bandwidth = 30;
            SNR = 20; % dB
        end
        
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
        
        function printSimulationProgress(outer_iter, inner_iter, max_outer_iter, max_inner_iter)
            ROW = max_outer_iter;
            COL = max_inner_iter;
            iRow = outer_iter;
            iCol = inner_iter;
            fprintf('Simulation progress\t%3.1f%%\n', ...
                100*(iCol+(iRow-1)*COL)/(ROW*COL) );
        end
		
	end

	
end
