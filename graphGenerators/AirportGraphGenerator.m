classdef AirportGraphGenerator < GraphGenerator
    
    properties % required by parent classes
		c_parsToPrint  = {'ch_name','s_edgeProbability','s_numberOfVertices'};
		c_stringToPrint  = {'','Edge prob.',''};
		c_patternToPrint = {'%s%s','%s = %g','%s%d vertices'};
	end

    
    properties(Constant)
        ch_name = 'Airport';
    end
    
    properties
        ch_filepath = '';   % path to data file
        relevantVariableName = {'ORIGIN_AIRPORT_ID', ...
            'DEST_AIRPORT_ID', 'DEP_DELAY', 'ARR_DELAY'};
    end
    
    methods
        % constructor
        function obj = AirportGraphGenerator(varargin)
            obj@GraphGenerator(varargin{:});
        end
        % realization: create graph
        function graph = realization(obj)
            T = obj.readDataFile();
            A = AirportGraphGenerator.createAdjacencyMatrix(T);
            graph = Graph('m_adjacency',A);
        end
    end
    
    methods 
    
        function T = readDataFile(obj)
            % READDATAFILE read airport data from csv/xlsx file
            %   Data file should be in the format specified by readtable function.
            % Since data can contain invalid items, NaN rows must be removed.
            filename = obj.ch_filepath;
            airportTable = readtable(filename);
            nanIndex = zeros(size(airportTable,1),1);
            for column = 1 : length(obj.relevantVariableName)
                nanIndex = nanIndex | ...
                    isnan(airportTable.(obj.relevantVariableName{column}));
            end
            airportTable(nanIndex,:) = [];  % delete relevant nan rows
            T = airportTable;
        end
    end
    
    methods(Static)
        
        function A = createAdjacencyMatrix(T)
            % create Laplacian matrix from table T
            % T is a table with each row denoting flight and delays
            % Node of graph is all the aiports
            % edges are # of flights between any two airports
            vertex = AirportGraphGenerator.getVertex(T); 
            A = zeros(length(vertex));

            ROWS = size(T,1);
            for r = 1 : ROWS
                origin = T.ORIGIN_AIRPORT_ID(r);
                dest = T.DEST_AIRPORT_ID(r);

                orgIndex = (vertex == origin);
                desIndex = (vertex == dest);
                A(orgIndex,desIndex) = A(orgIndex,desIndex) + 1;
            end
            A = A + A';
            A = double(A>0);
        end
        
       function VERTEX = getVertex(table)
            originAirport = unique(table.ORIGIN_AIRPORT_ID); % distinct
            destAirport = unique(table.DEST_AIRPORT_ID);
            VERTEX = union(originAirport, destAirport);
        end
    end
    
end

