function o = optimget2(options,name,default,flag)
%OPTIMGET2 Get OPTIM OPTIONS parameters for RSQP.
%   VAL = OPTIMGET2(OPTIONS,'NAME') extracts the value of the named parameter
%   from optimization options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  Case is ignored for 
%   parameter names.  [] is a valid OPTIONS argument.
%
%   VAL = OPTIMGET2(OPTIONS,'NAME',DEFAULT) extracts the named parameter as
%   above, but returns DEFAULT if the named parameter is not specified (is [])
%   in OPTIONS.
%
%   See also OPTIMGET, OPTIMSET2.

%   Copyright 2006 Institute of Industrial Control, Zhejiang University.
%   $Revision: 1.0 $  $Date: 2006/07/08 09:00:36 $

if nargin < 2
    error('MATLAB:optimget2:NotEnoughInputs', 'Not enough input arguments.');
end

try
    % get general optimization option parameters
	o = optimget(options,name,default,flag);
catch
    % get RSQP optimization option parameters
    if (nargin == 4) && isequal(flag,'fast')
        o = optimgetfast(options,name,default);
        return
    end
    names = {'basischange','basischoice','correction','lsmerit','preprocess','soc'};
    j = strmatch(lower(name),names)
    if isempty(j)
        error(sprintf(['Unrecognized property name ''%s''.  ' ...
                       'See OPTIMSET2 for possibilities.'], name))
    else
        o = options.(names{j});
        if isempty(o)
            o = default;
        end
    end
end

%--------------------------------------------------------------------------

function value = optimgetfast(options,name,defaultopt)

if ~isempty(options)
    value = options.(name);
else
    value = [];
end

if isempty(value)
    value = defaultopt.(name);
end

