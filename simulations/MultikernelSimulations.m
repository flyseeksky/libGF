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
			bl_estimator(1).c_replicatedVerticallyAlong = {'ch_name'};
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
			F = F_figure('X',Parameter.getXAxis(generator,sampler,est),'Y',mse,'leg',Parameter.getLegend(generator,sampler,est),'xlab',Parameter.getXLabel(generator,sampler,est));
			
		end
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%  simulations with MKL on synthetic data
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		
		% Simulation to test the regularization parameter
		function F = compute_fig_3001(obj,niter)
			
			SNR = 20; % dB
			N = 100;
						
			% generate graph
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability', 0.3,'s_numberOfVertices',N);
			graph = graphGenerator.realization();
			
			% generate Kernel matrix
			sigmaArray = linspace(0.01, 1.5, 20);
			nSigma = length(sigmaArray);
			m_kernel = zeros(N,N,nSigma);
			L = graph.getLaplacian();
			[V,D] = eig(L);
			d = diag(D);
			for iSigma = 1 : nSigma
				sigma = sigmaArray( iSigma );
				m_kernel(:,:,iSigma) = V * diag( exp( - sigma^2/2 * d ) ) * V';
			end
			uArray = logspace(-10, 0, 20);
			NMSE = nan(length(uArray));
			S = 50;
			
			% define graph function sampler
			sampler = UniformGraphFunctionSampler('s_numberOfSamples',S,'s_SNR',SNR);
			
			% geneartor graph function
			functionGenerator = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',30);
			m_graphFunction = functionGenerator.realization();
			
			[m_samples, m_positions] = sampler.sample(m_graphFunction);
			
			
			for iU = 1 : length(uArray)
				u = uArray(iU);
				
				estimator = MkrGraphFunctionEstimator('m_kernel', m_kernel, 's_mu', u);
				m_graphFunctionEstimate = estimator.estimate(m_samples, m_positions);
				NMSE(iU)  = norm(m_graphFunctionEstimate - m_graphFunction,'fro')^2/size(m_graphFunction,1);
				
				fprintf('Progress: %3.1f%%\n', ...
					100*( iU ) / ...
					( length(uArray)) );
			end
			
			% save NMSE.mat NMSE
			semilogx(uArray, NMSE)
			legend(sprintf('S = %d',S))
			
			
			F = [];
			
		end	
		
	end
	
	
end





