classdef TestLSEstimator < matlab.unittest.TestCase
    % Tests to be performed:
    %   1. illegal parameter: m_samples and m_positions not the same size
    %   2. illegal parameter: NaN values in paramter
    %   3. not engouh information: bandwidth not set
    %   4. not engouh information: graph not valid
    %   5. not consistent: m_samples and graph not the same size
    %   6. correct output: given graph, bandwidth, and parameter, valid
    %      corresponding output
    %   7. vector handling: given parameters, if bandwidth is set to be a
    %   vector, then perform estimation with each bandwith
    
    properties
        LSEstimator
        graph
        bandwidth = 5;
        N = 10;     % number of vertices
        p = 0.3;    % edge probability
        graphSignal
    end
    
    methods(TestMethodSetup)
        function createGraph(tc)
            ERGraphGenerator = ErdosRenyiGraphGenerator();
            ERGraphGenerator.s_numberOfVertices = tc.N;
            ERGraphGenerator.s_edgeProbability = tc.p;
            tc.graph = ERGraphGenerator.realization();
        end
        
        function createGraphSignal(tc)
            V = tc.graph.getLaplacianEigenvectors();
            alpha = zeros(tc.N,1);
            alpha(1:tc.bandwidth) = sort(rand(tc.bandwidth,1));
            tc.graphSignal = V*alpha;
        end
              
        function createLSEstimator(tc)
            tc.LSEstimator = LSGraphFunctionEstimator();
            tc.LSEstimator.graph = tc.graph;
            tc.LSEstimator.s_bandwidth = tc.bandwidth;
        end
    end
    
    methods(Test)
        function testGraphNotEmpty(tc)
            tc.verifyNotEmpty(tc.graph);
        end       
        function testGraphSize(tc)
            tc.verifySize(tc.graph.m_adjacency, [tc.N tc.N]);
        end
        function testGraphSignalBandwidth(tc)
            tc.assertNotEmpty(tc.graphSignal);
            V = tc.graph.getLaplacianEigenvectors();
            alpha = V'*tc.graphSignal;
            tc.verifyEqual(alpha(tc.bandwidth+1:end), ...
                zeros(tc.N-tc.bandwidth,1),'abstol',1e-5);
        end
        
        function testIllegalParameterSize(tc)
            m_positions = rand(tc.N,1);
            m_samples = rand(tc.N-1,1);
            tc.verifyError(@() tc.LSEstimator.estimate(m_samples, m_positions),...
                'LSEstimator:IllegalParameter');
        end
        
%         function testNanValues(tc)
%             S = 7;
%             m_positions = randperm(tc.N,S);
%             m_samples = tc.graphSignal(m_positions);
%             
%         end
    end
    
end

