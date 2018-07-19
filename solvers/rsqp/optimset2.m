function options = optimset2(varargin)
%OPTIMSET2 Create/alter optimization OPTIONS structure for RSQP.
%   OPTIONS = OPTIMSET2('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates an
%   optimization options structure OPTIONS in which the named parameters have
%   the specified values.  Any unspecified parameters are set to [] (parameters
%   with value [] indicate to use the default value for that parameter when
%   OPTIONS is passed to the optimization function).  Case is ignored for 
%   parameter names.
%
%OPTIMSET2 PARAMETERS for RSQP
%BasisChange - Basis change allowed during optimization [ {on} | off ]
%BasisChoice - Choose space decomposition strategy 
%              [ 'orthonormal' | 'orthogonal' | {'coordinate'} ]
%Correction - Calculate cross terms during optimization [ 0 | {1} | 2 ]
%LSMerit - Choose merit function for line search [ 0 | {1} ]
%preProcess - Dispose dependent linear equality constraints [ on | {off} ]
%SOC - Choose implementation strategy for the second order correction step
%      to overcome the Maratos effect [ 1 | {2} ]
%Other parameters see OPTIMSET for reference
%   See also OPTIMSET, OPTIMGET2.

%   Copyright 2006 Institute of Industrial Control, Zhejiang University.
%   $Revision: 1.0 $  $Date: 2006/07/07 22:02:25 $

% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
    optimset
    fprintf('            BasisChange: [ {on} | off ]\n');
    fprintf('            BasisChoice: [ orthonormal | orthogonal | {coordinate} ]\n');
    fprintf('             Correction: [ 0 | {1} | 2 ]\n');
    fprintf('                LSMerit: [ 0 | {1} ]\n');
    fprintf('             preProcess: [ on | {off} ]\n');
    fprintf('                    SOC: [ 1 | {2} ]\n');
    fprintf('\n');
    return;
end

% If pass in function name 'rsqp' then return the defaults.
if (nargin == 1) && strcmpi(varargin{1},'rsqp')
    options = optimset('fmincon');
    options.BasisChange = 'on';
    options.BasisChoice = 'coordinate';
    options.Correction = 1;
    options.LSMerit = 1;
    options.preProcess = 'off';
    options.SOC = 2;
    return
end

if mod(length(varargin),2)
    error('Arguments must occur in name-value pairs.');
end

rsqpopt = {'BasisChange','BasisChoice','Correction','LSMerit','preProcess','SOC'};
ind = [];
for i = 1:2:length(varargin)
    if mod(i,2) & find(strcmpi(varargin{i},rsqpopt))
        ind = [ind,i];
    end
end

% name-value pairs for general optimization options, varargin1; and
% name-value pairs for RSQP optimization options, varargin2
j = 1;
varargin1 = varargin;
for i = length(ind):-1:1
    varargin2{j} = varargin{ind(i)};
    varargin2{j+1} = varargin{ind(i)+1};
    j = j+2;
    varargin1(ind(i)+1) = [];
    varargin1(ind(i)) = [];    
end

options = optimset(varargin1{:});

if ~isempty(ind)
    for i = 1:2:length(varargin2)
        fname = varargin2{i};
        if ischar(varargin2{i+1})
            val = lower(varargin2{i+1});
        else
            val = varargin2{i+1};
        end
        [valid,errmsg] = checkfield(fname,val);
        if valid
            options.(fname) = val;
        else
            error(errmsg);
        end
    end
end

%-------------------------------------------------

function [valid, errmsg] = checkfield(field,value)

valid = 1;
errmsg = '';

% empty matrix is always valid
if isempty(value)
    return
end

switch field
    case {'BasisChange','preProcess'} % RSQP option: off, on
        if ~isa(value,'char') | ~any(strcmp(value,{'on';'off'}))
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
        end
    case {'BasisChoice'} % RSQP option: orthonormal, orthogonal, coordinate
        if ~isa(value,'char') | ~any(strcmp(value,{'orthonormal';'orthogonal';'coordinate'}))
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''orthonormal'', ''orthogonal'' or ''coordinate''.',field);
        end
    case {'Correction'} % RSQP option: 0, 1, 2
        if ~isa(value,'double') | ~isequal(value,0) & ~isequal(value,1) & ~isequal(value,2)
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be 0, 1 or 2.',field);
        end
    case {'LSMerit'} % RSQP option: 0, 1
        if ~isa(value,'double') | ~isequal(value,0) & ~isequal(value,1)
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be 0 or 1.',field);
        end
    case {'SOC'} % RSQP option: 1, 2
        if ~isa(value,'double') | ~isequal(value,1) & ~isequal(value,2)
            valid = 0;
            errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be 1 or 2.',field);
        end
    otherwise
        valid = 0;
        error('Unknown field name for Options structure.')
end

