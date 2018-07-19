function result = proxy(fname,outnum,varargin)
% interface of RSQP Toolbox and Optimization Toolbox

%   Copyright 2006 Institute of Industrial Control, Zhejiang University.
%   $Revision: 1.0 $  $Date: 2006/07/12 19:25:49 $

result = cell(outnum,1);
try
	[result{:}] = feval(fname,varargin{:});
catch
	error(lasterr)
end

