%
%  FIGURES FOR THE PAPER ON MULTIKERNEL
%
% 

classdef MultikernelSimulations < simGFFunctionSet
	
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
		
	end
	
	
end





