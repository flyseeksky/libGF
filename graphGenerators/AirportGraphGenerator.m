classdef AirportGraphGenerator < GraphGenerator
    
    properties % required by parent classes
        c_parsToPrint  = {'ch_name','s_numberOfVertices'};
        c_stringToPrint  = {'',''};
        c_patternToPrint = {'%s%s','%s%d vertices'};
    end
    
    
    properties(Constant)
        ch_name = 'Airport';
    end
    
    properties
        ch_filepath = '';   % path to data file
        airportTable;       % table that contains raw data
        delayTable;         % a table that each entry contains the average 
                            % delay of one day and one airport
                            
        v_listOfAirport;    % valid airports ID
        v_listOfDay;        % days involved
        
        c_relevantVariableName = {'ORIGIN_AIRPORT_ID', ...
            'DEST_AIRPORT_ID', 'DEP_DELAY', 'ARR_DELAY'};   % variables
                            % that we care about
    end
    
    methods
        % constructor
        function obj = AirportGraphGenerator(varargin)
            obj@GraphGenerator(varargin{:});
        end
        
        % realization: create graph
        function graph = realization(obj)         
                % read data into a table
                obj.airportTable = AirportGraphGenerator.importAirportData(obj.ch_filepath);
                
                %% remove all NaN entries
                obj = obj.removeNanRows();
                
                %% average all fights in the same day and the same date
                [depDelayMatrix, arrDelayMatrix] = obj.getDelayTable();
                %% create adjacency matrix
                
                % create graph
                adjacency = AirportGraphGenerator.createAdjacencyMatrix(obj.airportTable);
                graph = Graph('m_adjacency',adjacency);
        end
    end
    
    methods  % getter and setter
        function listOfAirport = get.v_listOfAirport(obj)
            if isempty(obj.v_listOfAirport)
                airportIDList = obj.airportTable.ORIGIN_AIRPORT_ID;
                listOfAirport = sort(unique(airportIDList));  % ascend
            else
                listOfAirport = obj.v_listOfAirport;
            end
        end
        
        function listOfDay = get.v_listOfDay(obj)
            if isempty(obj.v_listOfDay)
                dayList = obj.airportTable.FL_DATE;
                listOfDay = sort(unique(dayList));            % ascend
            else
                listOfDay = obj.v_listOfDay;
            end
        end
    end
    
    methods
        function [depDelayMatrix, arrDelayMatrix] = getDelayTable(obj)
            % get the two tables from raw data
            % depDelayTable:
            %       each entry contains the average departure delay of 
            %       that day and that airport
            % depDelayTable:
            %       each entry contains the average arrival delay of 
            %       that day and that airport
            listOfDay = obj.v_listOfDay;
            listOfAirport = obj.v_listOfAirport;
            
            depDelayMatrix = NaN(length(listOfAirport), length(listOfDay));
            arrDelayMatrix = NaN(length(listOfAirport), length(listOfDay));
            
            for day = 1:length(listOfDay)
                % create a sub-table that contains rows of that day
                dayTable = obj.airportTable( obj.airportTable.FL_DATE == listOfDay(day), : );
                for airport = 1:length(listOfAirport)
                    % create a sub-table that contains rows of that airport
                    airportDayTable = dayTable( dayTable.ORIGIN_AIRPORT_ID == ...
                        listOfAirport(airport) ,:);

                        depDelayMatrix(airport, day) = mean(airportDayTable.DEP_DELAY);
                        arrDelayMatrix(airport, day) = mean(airportDayTable.ARR_DELAY);
                end
            end
        end
        
        function obj = removeNanRows(obj)
            % remove all nan rows that contains in obj.airportTable
            rows = false( size(obj.airportTable,1),1 );
            for col = 1:length(obj.c_relevantVariableName)
                varCol = obj.airportTable.(obj.c_relevantVariableName{col});
                rows = rows | isnan(varCol);
            end
            obj.airportTable(rows, :) = [];     % delete nan rows
        end
        
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
        
        function March2016 = importAirportData(filename, startRow, endRow)
            %IMPORTFILE Import numeric data from a text file as a matrix.
            %   MARCH2016 = IMPORTAIRPORTDATA(FILENAME) Reads data from text file FILENAME for
            %   the default selection.
            %
            %   MARCH2016 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
            %   STARTROW through ENDROW of text file FILENAME.
            %
            % Example:
            %   March2016 = importfile('March2016.csv', 2, 479123);
            %
            %    See also TEXTSCAN.
            
            % Auto-generated by MATLAB on 2016/05/19 18:28:31
            
            %% Initialize variables.
            delimiter = ',';
            if nargin<=2
                startRow = 2;
                endRow = inf;
            end
            
            %% Format string for each line of text:
            %   column1: datetimes (%{MM/dd/yyyy}D)
            %	column2: double (%f)
            %   column3: double (%f)
            %	column4: double (%f)
            %   column5: text (%s)
            %	column6: double (%f)
            %   column7: text (%s)
            %	column8: double (%f)
            %   column9: text (%s)
            %	column10: double (%f)
            %   column11: text (%s)
            %	column12: text (%s)
            %   column13: text (%s)
            %	column14: text (%s)
            %   column15: text (%s)
            %	column16: text (%s)
            % For more information, see the TEXTSCAN documentation.
            formatSpec = '%{MM/dd/yyyy}D%f%f%f%s%f%s%f%s%f%s%s%s%s%s%s%[^\n\r]';
            
            %% Open the text file.
            fileID = fopen(filename,'r');
            
            %% Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col};dataArrayBlock{col}];
                end
            end
            
            %% Close the text file.
            fclose(fileID);
            
            %% Post processing for unimportable data.
            % No unimportable data rules were applied during the import, so no post
            % processing code is included. To generate code which works for
            % unimportable data, select unimportable cells in a file and regenerate the
            % script.
            
            %% Create output variable
            March2016 = table(dataArray{1:end-1}, 'VariableNames', ...
                {'FL_DATE','AIRLINE_ID','FL_NUM','ORIGIN_AIRPORT_ID',...
                'ORIGIN','DEST_AIRPORT_ID','DEST','DEP_DELAY',...
                'DEP_TIME_BLK','ARR_DELAY','ARR_TIME_BLK',...
                'CARRIER_DELAY','WEATHER_DELAY','NAS_DELAY',...
                'SECURITY_DELAY','LATE_AIRCRAFT_DELAY'});
            
            % For code requiring serial dates (datenum) instead of datetime, uncomment
            % the following line(s) below to return the imported dates as datenum(s).
            
            % March2016.FL_DATE=datenum(March2016.FL_DATE);
            
            
        end
        
    end
    
end