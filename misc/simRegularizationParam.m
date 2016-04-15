function simRegularizationParam

MONTE_CARLO = 200;
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

end