classdef RealDataGraphFunctionGenerator  < GraphFunctionGenerator
	
	
	properties % Required by superclass Parameter
		c_parsToPrint    = {};
		c_stringToPrint  = {};
		c_patternToPrint = {};
	end 
	
	properties
		name = 'RealSignal';
		v_realSignal% integer		
	end
	
	methods
		
		function obj = RealDataGraphFunctionGenerator(varargin)
			% constructor
			obj@GraphFunctionGenerator(varargin{:});
		end
		
		
		function M_graphFunction = realization(obj,s_numberOfRealizations)
			% M_GRAPHFUNCTION   N x S_NUMBEROFREALIZATIONS matrix where N is
			%                   the number of vertices.The signal is copied
			%                   s_numberOfRealizations times
			
			assert(~isempty(obj.graph));

			if nargin < 2
				s_numberOfRealizations = 1;
			end
			
            M_graphFunction=repmat(obj.v_realSignal,[1,s_numberOfRealizations]);
			v_powerVect=sqrt(sum(M_graphFunction.^2,1));
			M_graphFunction=M_graphFunction./repmat(v_powerVect,size(M_graphFunction,1),1);
			
        end
		
		
	end
	
end

