classdef ExponentiallyDecayingFunctioinGenerator < GraphFunctionGenerator
    % Generate exponentially decaying function on graphs
    % need to specify the graph and bandwidth
    
    properties
        graph
        s_bandwidth
    end
    
    methods
        function obj = ExponentiallyDecayingFunctioinGenerator(vargin)
            obj@GraphFunctionGenerator(vargin);
        end
        
        function M_graphFunction = realization(obj,s_numberOfRealizations)
            if isempty(obj.graph) || isempty(obj.s_bandwidth)
                error('ExponentiallyDecayingFunctionGenerator: Paramter not set correctly')
            end
            
            L = graph.getLaplacian();
        end
    end
    
end

