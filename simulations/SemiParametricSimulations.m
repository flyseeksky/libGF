%
%  FIGURES FOR THE PAPER ON MULTIKERNEL
%
%

classdef SemiParametricSimulations < simFunctionSet
	
	properties
		
	end
	
	methods
		
		
		function F = compute_fig_1001(obj,niter)
			
			%0. define parameters
			s_sigma=1.3;
			s_numberOfClusters=5;
			s_lambda=10^-5;
			s_monteCarloSimulations=niter;%100;
			s_bandwidth1=10;
			s_bandwidth2=20;
			s_SNR=1000;
			s_beta=0.02;
			s_alpha=0.005;
			s_niter=10;
			s_epsilon=0.2;
			s_functionTypes=5;
			v_sampleSetSize=(0.1:0.1:1);
			
			m_meanSquaredError=zeros(size(v_sampleSetSize,2),s_functionTypes);
			% define graph
			[Ho,Mo,Alto,Hn,Mn,Altn] = readTemperatureDataset;
			%tic
			graphGenerator = GraphLearningSmoothSignalGraphGenerator('m_observed',Ho,'s_niter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta,'m_missingValuesIndicator',[]);
			%toc
			
			% tic
			% graphGenerator=SmoothSignalGraphGenerator('m_observed',Ho,'s_maxIter',s_niter,'s_alpha',s_alpha,'s_beta',s_beta);
			% toc
			
			graph = graphGenerator.realization;
			%graph1 = graphGenerator1.realization;
			%L1=graph.getLaplacian
			%L2=graph1.getLaplacian
			v_sampleSetSize=round(v_sampleSetSize*graph.getNumberOfVertices);
			
			
			m_basis= SemiParametricSimulations.parametricPartForTempData(graph.getLaplacian,Alto,s_numberOfClusters);
			
			%functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth);
			%functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
			%signal
			v_realSignal=Hn(:,3);
			functionGenerator=RealDataGraphFunctionGenerator('graph',graph,'v_realSignal',v_realSignal);
			% define bandlimited function estimator
			m_laplacianEigenvectors=(graph.getLaplacianEigenvectors);
			bandlimitedGraphFunctionEstimator1 = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors',m_laplacianEigenvectors(:,1:s_bandwidth1));
			bandlimitedGraphFunctionEstimator2 = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors',m_laplacianEigenvectors(:,1:(s_bandwidth2)));
			
			% define Kernel function
			diffusionGraphKernel = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix);
			
			%define semi-parametric estimator
			semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			
			% Simulation
			for s_sampleSetIndex=1:size(v_sampleSetSize,2)
				
				
				s_numberOfSamples=v_sampleSetSize(s_sampleSetIndex);
				%sample
				sampler = UniformGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR);
				m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
				[m_samples,m_positions] = sampler.sample(m_graphFunction);
				%estimate
				m_graphFunctionEstimateBL1 = bandlimitedGraphFunctionEstimator1.estimate(m_samples,m_positions);
				m_graphFunctionEstimateBL2= bandlimitedGraphFunctionEstimator2.estimate(m_samples,m_positions);
				m_graphFunctionEstimateNP=nonParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEstimateSP=semiParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,s_epsilon);
				
				% Performance assessment
				m_meanSquaredError(s_sampleSetIndex,1) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateBL1,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,2) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateBL2,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,3) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateNP,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,4) = SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateSP,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,5) = SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEpsilonInsesitiveEstimateSP,m_graphFunction);
				
			end
			%m_meanSquaredError(m_meanSquaredError>1)=1;
			%save('real.mat');
			
			F = F_figure('X',v_sampleSetSize,'Y',m_meanSquaredError','xlab','Number of observed vertices (S)','ylab','NMSE','leg',{strcat('Bandlimited  ',sprintf(' W=%g',s_bandwidth1)),strcat('Bandlimited',sprintf(' W=%g',s_bandwidth2)),'Nonparametric (SL)','Semi-parametric (SL)','Semi-parametric (\epsilon-IL)'});
		end
		function F = compute_fig_1002(obj,niter)
			F = obj.load_F_structure(1001);
			F.ylimit=[0 1];
			F.xlimit=[10 89];
			F.styles = {'-','--','-o','-x','--^'};
			F.pos=[680 729 509 249];
			F.leg={strcat('Bandlimited  ',sprintf(' W=10')),strcat('Bandlimited',sprintf(' W=20')),'Nonparametric (SL)','Semi-parametric (SL)','Semi-parametric (\epsilon-IL)'};

			F.leg_pos = 'north';      % it can be 'northwest',
			%F.leg_pos_vec = [0.547 0.673 0.182 0.114];
		end
		function F = compute_fig_1003(obj,niter)
			
			%Assilomar Synthetic simmulation
			
			%0. define parameters
			s_sigma=0.7;
			s_numberOfClusters=4;
			s_Type=2;
			s_lambda=10^-8;
			s_monteCarloSimulations=niter;
			s_bandwidth1=10;
			s_bandwidth2=20;
			s_SNR=4;
			s_dataSetSize=100;
			s_functionTypes=5;
			s_epsilon=0.1;
			v_sampleSetSize=round((0.1:0.1:1)*s_dataSetSize);
			m_meanSquaredError=zeros(size(v_sampleSetSize,2),s_functionTypes);
			% define graph function generator Parametric basis
			graphGenerator = ErdosRenyiGraphGenerator('s_edgeProbability',.6,'s_numberOfVertices',s_dataSetSize);
			graph = graphGenerator.realization;
			
			m_sparseBasis=graph.getClusters(s_numberOfClusters,s_Type);
			m_basis=1.5*full(m_sparseBasis);
			%m_basis=m_basis*1;
			
			functionGeneratorBL = BandlimitedGraphFunctionGenerator('graph',graph,'s_bandwidth',s_bandwidth1);
			functionGenerator= SemiParametricGraphFunctionGenerator('graph',graph,'graphFunctionGenerator',functionGeneratorBL,'m_parametricBasis',m_basis);
			m_laplacianEigenvectors=(graph.getLaplacianEigenvectors);
			
			% define bandlimited function estimator
			bandlimitedGraphFunctionEstimator1 = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors',functionGeneratorBL.basis);
			bandlimitedGraphFunctionEstimator2 = BandlimitedGraphFunctionEstimator('m_laplacianEigenvectors',m_laplacianEigenvectors(:,1:(s_bandwidth2)));
			% define Kernel function
			diffusionGraphKernel = DiffusionGraphKernel('m_laplacian',graph.getLaplacian,'s_sigma',s_sigma);
			%define non-parametric estimator
			nonParametricGraphFunctionEstimator=NonParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix);
			
			%define semi-parametric estimator
			semiParametricGraphFunctionEstimator = SemiParametricGraphFunctionEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			%define semi-parametric epsinlon  insensitive estimator
			semiParametricGraphFunctionEpsilonInsesitiveEstimator=SemiParametricGraphFunctionEpsilonInsesitiveEstimator('m_kernels',diffusionGraphKernel.generateKernelMatrix,'m_basis',m_basis);
			
			% Simulation
			for s_sampleSetIndex=1:size(v_sampleSetSize,2)
				
				
				s_numberOfSamples=v_sampleSetSize(s_sampleSetIndex);
				%sample
				sampler = UniformGraphFunctionSampler('s_numberOfSamples',s_numberOfSamples,'s_SNR',s_SNR);
				m_graphFunction = functionGenerator.realization(s_monteCarloSimulations);
				[m_samples,m_positions] = sampler.sample(m_graphFunction);
				%estimate
				m_graphFunctionEstimateBL1 = bandlimitedGraphFunctionEstimator1.estimate(m_samples,m_positions);
				m_graphFunctionEstimateBL2= bandlimitedGraphFunctionEstimator2.estimate(m_samples,m_positions);
				m_graphFunctionEstimateNP=nonParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEstimateSP=semiParametricGraphFunctionEstimator.estimate(m_samples,m_positions,s_lambda);
				m_graphFunctionEpsilonInsesitiveEstimateSP=semiParametricGraphFunctionEpsilonInsesitiveEstimator.estimate(m_samples,m_positions,s_lambda,s_epsilon);
				% Performance assessment
				m_meanSquaredError(s_sampleSetIndex,1) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateBL1,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,2) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateBL2,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,3) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateNP,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,4) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEstimateSP,m_graphFunction);
				m_meanSquaredError(s_sampleSetIndex,5) =SemiParametricSimulations.estimateMeanSquaredError(m_graphFunctionEpsilonInsesitiveEstimateSP,m_graphFunction);
			end
			%m_meanSquaredError(m_meanSquaredError>1)=1;
			%save('synthetic.mat');
			F = F_figure('X',v_sampleSetSize,'Y',m_meanSquaredError','xlab','Number of observed vertices (S)','ylab','NMSE','leg',{strcat('Bandlimited  ',sprintf(' W=%g',s_bandwidth1)),strcat('Bandlimited',sprintf(' W=%g',s_bandwidth2)),'Nonparametric (SL)','Semiparametric (SL)','Semiparametric (\epsilon-IL)'});
		
		
		end
		function F = compute_fig_1004(obj,niter)
			F = obj.load_F_structure(1003);
			F.ylimit=[0 1];
			F.xlimit=[10 100];
			F.styles = {'-','--','-o','-x','--^'};
			F.pos=[680 729 509 249];
			F.leg={strcat('Bandlimited  ',sprintf(' W=10')),strcat('Bandlimited',sprintf(' W=20')),'Nonparametric (SL)','Semi-parametric (SL)','Semi-parametric (\epsilon-IL)'};

			%F.leg_pos = 'northeast';      % it can be 'northwest',
			F.leg_pos_vec = [0.647 0.683 0.182 0.114];
		end
	end
	
	
	
	methods(Static)
		
		function B=parametricPartForTempData(L,feat,n_clusters)
			
			
			C=kmeans(feat,n_clusters);
			% cluster via the altitude information..
			B=zeros(size(L,1),1);
			for i=1:n_clusters
				B(C==i,i)=1;
			end
			
		end
		
		
		
		function res=estimateMeanSquaredError(m_est,m_observed)
			res=0;
			for i=1:size(m_est,2)
				res=res+norm(m_est(:,i)-m_observed(:,i))^2/norm(m_observed(:,i))^2;
			end
			res=(1/size(m_est,2))*res;
		end
		
	end
	
	
end





