function NMSE = Simulate( generator, sampler, ...
    estimator, MONTE_CARLO, refSignal)
% SIMULATE      simulate function estiamtion on graph
% Input:
%       generator       graph function generator matrix
%       sampler         graph function sampler matrix
%       estimator       graph function estiamtor matrix
%       MONTE_CARLO     number of Monte Carlo simulations
% Output:
%       NMSE            normalized mean squared error, whose size is
%                       consistent with generator/sampler/estimator
%

DEBUG = true;

% get the dimension of result
ROW = max( [size(generator,1), size(sampler,1), size(estimator,1)] );
COL = max( [size(generator,2), size(sampler,2), size(estimator,2)] );

% simulation
NMSE = NaN(ROW,COL);
normRef = norm(refSignal);
for iRow = 1 : ROW
    for iCol = 1 : COL
        m_estimate = MonteCarloSimulation( ...
            getObjectFromMat(generator,iRow, iCol), ...
            getObjectFromMat(sampler,iRow,iCol), ...
            getObjectFromMat(estimator, iRow, iCol), ...
            MONTE_CARLO);
        NMSE(iRow, iCol) = norm(m_estimate - refSignal) / normRef;
        
        if DEBUG
            fprintf('Simulation progress\t%3.1f%%\n', ...
                100*(iCol+(iRow-1)*COL)/(ROW*COL) );
        end
    end
end

end

function obj = getObjectFromMat(objMat, row, col)
% this function is used to retrieve the obj (generator/sampler/estiamtor)
% from replicated matrix by specifiying row and col.
% If there is only 
[M,N] = size(objMat);

if M == 1 && N == 1
    % only one object in objMat
    % always return this object
    obj = objMat;
elseif M == 1
    % objMat is a row vector, only column index matters
    obj = objMat(col);
elseif N == 1
    % objMat is a column vector, only row index matters
    obj = objMat(row);
else
    % objMat is a non-trivial matrix, return corresponding element
    obj = objMat(row, col);
end

end

function m_estimate = MonteCarloSimulation( generator, sampler, estimator, MONTE_CARLO )

graphFunction = generator.realization();
parfor iMonte = 1 : MONTE_CARLO
    [m_samples, m_positions] = sampler.sample(graphFunction);
    estimate(:, iMonte) = estimator.estimate(m_samples, m_positions);
end

m_estimate = mean(estimate,2);

end