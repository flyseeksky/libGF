classdef Parameter 
	% This class allows obtaining information about parameters of
	% descendants (e.g. generators, samplers, estimators) in text form.
	% This is useful to make figures. 
	%
	% all classes extending Parameter shall allow their constructors to be
	% called without any parameter.
        %
        % 

	properties(Abstract)
		% the n-th property of an object of class Parameter will be printed
		% as   my_sprintf( obj.c_patternToPrint{n} , obj.c_stringToPrint{n} , getfield(obj,obj.c_parsToPrint{n}) );		
		c_parsToPrint    % cell array of strings containing the name of the parameters to print
		c_stringToPrint  % cell array of strings containing the string for printing each parameter (titles, axes and legends)
		c_patternToPrint % cell array of strings containing the string pattern for sprintf in titles (e.g. {'%s = %d'})
	end 
	
	properties
		replicated_vertically_along = {};
		replicated_horizontally_along = {};
	end
	
	properties(Constant)
		def_chars_per_line = 200;
	end
	
	methods
		
		function obj = Parameter(varargin)
			% varargin is a sequence of (property_name,property_value)
			% pairs
			% Example:
			%     car = Parameter('numberOfWheels',4,'year',2009,'color','blue')
			obj =  assignParametersByName(obj,varargin{:});
			
		end
		
		
		function str = getParameterByNumber( obj1, par_number )
			obj1 = obj1(1,1);
			assert(numel(obj1)==1);			
			str = my_sprintf( obj1.c_patternToPrint{par_number} , obj1.c_stringToPrint{par_number} , my_getfield(obj1,obj1.c_parsToPrint{par_number}) );			
        end
		function str = getParameterByName( obj1, name )
			par_number = obj1.parameterIndex(name);
			str =  getParameterByNumber( obj1, par_number );
		end
		

		function ind = parameterIndex(obj,str)
			for k=1:length(obj.c_parsToPrint)
				if strcmp(obj.c_parsToPrint{k},str)
					ind = k;
					return
				end
			end
			str = sprintf('parameter ''%s'' does not exist in %s.c_parsToPrint',str,class(obj));
			error(str);
		end
		
	end

	
	methods
				
		function obj_mat = replicate(obj,fieldname_1,fieldvalues_1,fieldname_2,fieldvalues_2)
			% input 
			%    FIELDVALUES_1: cell array with M values for the file in
			%                   Parameter called FIELDNAME_1
			%    FIELDNAME_1  : string with the name of the field
			%    
			%    FIELDVALUES_2: cell array with N values for the file in
			%                   Parameter called FIELDNAME_1
			%    FIELDNAME_2  : string with the name of the field
			%
			% output
			%    obj_mat      : MxN matrix where the (m,n)-th element has
			%                   all the values equal to those of the
			%                   current object except from the field
			%                   FIELDNAME_1, which has value
			%                   FIELDVALUES_1{m} and the field
			%                   FIELDNAME_2, which has value
			%                   FIELDVALUES_2{n}, that is, FIELDNAME_1
			%                   replicates the object vertically, and
			%                   FIELDNAME_2 does it horizontally
			%
			
			assert(numel(obj)==1,'REPLICATE not implemented for input array objects');
			% this  function is implemented for scalar objects. see
			% replicate_horizontally above
			
			% check that fieldvalues are cells. You can use num2cell
			assert( (isempty(fieldvalues_1)||iscell(fieldvalues_1))&&(isempty(fieldvalues_2)||iscell(fieldvalues_2) ));
			
			obj.replicated_vertically_along={fieldname_1};
			obj.replicated_horizontally_along={fieldname_2};
			M = max(length(fieldvalues_1),1);
			N = max(length(fieldvalues_2),1);
			obj_mat = repmat(obj,M,N);
			for m=1:M
				for n=1:N
					obj_mat(m,n)=obj;
					if ~isempty(fieldname_1)
						obj_mat(m,n).(fieldname_1) = fieldvalues_1{m};
					end
					if ~isempty(fieldname_2)
						obj_mat(m,n).(fieldname_2) = fieldvalues_2{n};
					end
				end
			end
			
		end
				
		function d = is_replicated_vertically_along(obj,str)
			
			obj = obj(1,1);
			
			if isempty( obj.replicated_vertically_along )
				d = 0;
			else
				d = 0;
				for k=1:length( obj.replicated_vertically_along )
					if strcmp( str , obj.replicated_vertically_along{k} )
						d=1;
						return
					end
				end
			end
			
		end
				
		function d = is_replicated_horizontally_along(obj,str)
			
			obj = obj(1,1);
			
			if isempty( obj.replicated_horizontally_along )
				d = 0;
			else
				d = 0;
				for k=1:length( obj.replicated_horizontally_along )
					if strcmp( str , obj.replicated_horizontally_along{k} )
						d=1;
						return
					end
				end
			end
			
		end
				
	end
	
	
	methods(Static)
		

		
		
		function leg = formLegend(varargin)
			
			list = Parameter.formLegendList(varargin{:});
			
			leg ={};
			for k = 1:size(list,1)
				leg{k} = ParameterArray.strListToText( {list{k,:}},1000 );
			end
			
		end
		
		function leg = formLegendList(varargin)
			% it makes a 2D cell array with the entries for the legend
			col = 1;
			leg = {};
			for k=1:nargin
				obj_now = varargin{k};
				obj11 = obj_now(1,1);
				leg_pars = obj11.replicated_vertically_along;
				
				for par = 1:length(leg_pars)
					par_now = leg_pars{par};
					if isempty(par_now)
						continue
					end
					for row = 1:size(obj_now,1)
						str = obj_now(row,1).getParameterByName(par_now);
						leg{row,col} = str;
					end
					
					col = col + 1;
				end
			end
		end
				
		function tit = formTitleFilter(no_list,varargin)
			global chars_per_line
			if isempty(chars_per_line)
				chars_per_line = Parameter.def_chars_per_line;
			end
			
			list = Parameter.formTitleList(varargin{:});
			assert(max(no_list)<=length(list));
			inds = setxor(1:length(list),no_list);
			list = list(inds);
			tit = [];
			tit_len = 0;
			for k=1:length(list)
				% end of line
				if tit_len + length(list{k}) > chars_per_line
					tit = [tit sprintf('\n')];
					tit_len = length(list{k});										
				else					
					tit_len = tit_len + length(list{k});					
				end
				tit = [tit list{k}];
				% comma
				if k~=length(list)
					tit = [tit ', '];
				end
			end

		end
				
		function list = formTitleList(varargin)
			% cell array with the strings to put in the title
			% arguments are objects of the class Parameter

			list = {};
			for nn = 1:nargin
				obj1 = varargin{nn};				
				obj1 = obj1(1,1);
				for k=1:length(obj1.c_parsToPrint)
					par_name = obj1.c_parsToPrint{k};
					if (~is_replicated_vertically_along(obj1,par_name) )&&(~is_replicated_horizontally_along(obj1,par_name) )
						str = obj1.getParameterByNumber(k);
						if ~isempty(str)
							list = {list{:},str};
						end
					end
				end
			end
		end
				
		function tit = formTitle(varargin)
			global chars_per_line
			if isempty(chars_per_line)
				chars_per_line = Parameter.def_chars_per_line;
			end
			
			list = Parameter.formTitleList(varargin{:});
			tit = ParameterArray.strListToText(list,chars_per_line);

		end
				
	end
	
	
end

