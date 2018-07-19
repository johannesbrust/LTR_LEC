function [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = rsqp(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
%RSQP finds a constrained minimum of a function of large number of variables.
%   RSQP solves problems of the form:
%       min F(X)  subject to:  Aeq*X  = Beq (linear constraints)
%        X                     Ceq(X) = 0  (nonlinear constraints)
%                              LB <= X <= UB                       
%
%   X=RSQP(FUN,X0,Aeq,Beq) minimizes FUN subject to the linear equalities
%   Aeq*X = Beq.
%
%   X=RSQP(FUN,X0,Aeq,Beq,LB,UB) defines a set of lower and upper
%   bounds on the design variables, X, so that a solution is found in 
%   the range LB <= X <= UB. Use empty matrices for LB and UB
%   if no bounds exist. Set LB(i) = -Inf if X(i) is unbounded below; 
%   set UB(i) = Inf if X(i) is unbounded above.
%
%   X=RSQP(FUN,X0,Aeq,Beq,LB,UB,NONLCON) subjects the minimization to the 
%   constraints defined in NONLCON. The function NONLCON accepts X and returns 
%   the vectors Ceq, representing the nonlinear equalities. RSQP minimizes FUN
%   such that Ceq(X)=0. 
%   (Set LB=[] and/or UB=[] if no bounds exist.)
%
%   X=RSQP(FUN,X0,Aeq,Beq,LB,UB,NONLCON,OPTIONS) minimizes with the 
%   default optimization parameters replaced by values in the structure OPTIONS, 
%   an argument created with the OPTIMSET2 function.  See OPTIMSET2 for details.  Used
%   options are Display, TolX, TolFun, TolCon, DerivativeCheck, Diagnostics, GradObj, 
%   GradConstr, Hessian, MaxFunEvals, MaxIter, DiffMinChange and DiffMaxChange, 
%   LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX, Hessian, HessMult, 
%   HessPattern, BasisChange, BasisChoice, Correction, LSMerit, preProcess, SOC. 
%   Use the GradObj option to specify that FUN also returns a second output argument G 
%   that is the partial derivatives of the function df/dX, at the point X. 
%   Use the Hessian option to specify that FUN also returns a third output argument H 
%   that is the 2nd partial derivatives of the function (the Hessian) at the point X.  
%   The Hessian is not used by the line-search method. Use the GradConstr 
%   option to specify that NONLCON also returns third and fourth output arguments 
%   GC and GCeq, where GC is the partial derivatives of the constraint vector of 
%   inequalities C, and GCeq is the partial derivatives of the constraint vector of 
%   equalities Ceq (NOTE: there is no nonlinear inequalities C in RSQP).
%   Use the BasisChange option to specify whether basis change is allowed during 
%   optimization. Use the BasisChoice option to specify space decomposition strategy. 
%   Use the Correction option to specify whether to calculate the cross terms. 
%   Use the LSMerit option to specify merit function for line search.
%   Use the preProcess option to specify whether to treat with dependent linear
%   equality constraints at the beginning of optimization.
%   Use the SOC option to specify implementation strategy for the second order 
%   correction step to overcome the Maratos effect
%   Use OPTIONS = [] as a place holder if no options are set.
%  
%   X=RSQP(FUN,X0,Aeq,Beq,LB,UB,NONLCON,OPTIONS,P1,P2,...) passes the 
%   problem-dependent parameters P1,P2,... directly to the functions FUN 
%   and NONLCON: feval(FUN,X,P1,P2,...) and feval(NONLCON,X,P1,P2,...).  Pass   
%   empty matrices for Aeq, Beq, OPTIONS, LB, UB, and NONLCON to use the 
%   default values.
%
%   [X,FVAL]=RSQP(FUN,X0,...) returns the value of the objective 
%   function FUN at the solution X.
%
%   [X,FVAL,EXITFLAG]=RSQP(FUN,X0,...) returns a string EXITFLAG that 
%   describes the exit condition of RSQP. Possible values of EXITFLAG and the corresponding 
%   exit conditions are
%
%     1  First order optimality conditions satisfied to the specified tolerance.
%     2  Change in X less than the specified tolerance.
%     3  Change in the objective function value less than the specified tolerance.
%     4  Magnitude of search direction smaller than the specified tolerance and 
%         constraint violation less than options.TolCon.
%     5  Magnitude of directional derivative less than the specified tolerance
%         and constraint violation less than options.TolCon.
%     0  Maximum number of function evaluations or iterations reached.
%    -1  Optimization terminated prematurely by user.
%    -2  No feasible point found.
%    -3  CPU time limit was exceeded
%    -4  NaN occurred
%   
%   [X,FVAL,EXITFLAG,OUTPUT]=RSQP(FUN,X0,...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the number
%   of function evaluations in OUTPUT.funcCount, the algorithm used in 
%   OUTPUT.algorithm, the number of CG iterations (if used) in OUTPUT.cgiterations, 
%   and the first-order optimality (if used) in OUTPUT.firstorderopt.
%
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA]=RSQP(FUN,X0,...) returns the Lagrange multipliers
%   at the solution X: LAMBDA.lower for LB, LAMBDA.upper for UB, LAMBDA.eqlin is for 
%   the linear equalities,and LAMBDA.eqnonlin is for the nonlinear equalities.
%   
%   [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD]=RSQP(FUN,X0,...) returns the value of 
%   the gradient of FUN at the solution X.
%
%   Examples
%       FUN can be specified using @:
%           X = rsqp(@humps,...)
%       In this case, F = humps(X) returns the scalar function value F of the HUMPS function
%       evaluated at X.
%
%       FUN can also be an inline object:
%           X = rsqp(inline('3*sin(x(1))+exp(x(2))'),[1;1],[],[],[0 0])
%           returns X = [0;0].
%
%   See also OPTIMSET2, FMINCON, FMINUNC, FMINBND, FMINSEARCH, @, INLINE.

%   Copyright 2006 Institute of Industrial Control, Zhejiang University.
%   $Revision: 1.0 $  $Date: 2006/07/08 15:23:12 $

defaultopt = struct('Display','iter','LargeScale','on', ...
    'TolX',1e-6,'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off',...
    'Diagnostics','off','FunValCheck','off',...
    'GradObj','off','GradConstr','off',...
    'HessMult',[],...   % HessMult [] by default
    'Hessian','off','HessPattern','sparse(ones(numberOfVariables))',...
    'MaxFunEvals','100*numberOfVariables',...
    'MaxSQPIter',inf,...
    'DiffMaxChange',1e-1,'DiffMinChange',1e-8,...
    'PrecondBandWidth',0,'TypicalX','ones(numberOfVariables,1)',...
    'MaxPCGIter','max(1,floor(numberOfVariables/2))',...
    'TolPCG',0.1,'MaxIter',400,...
    'BasisChange','on',...              % option for RSQP
    'BasisChoice','coordinate',...      % option for RSQP
    'Correction',1,...                  % option for RSQP
    'LSMerit',1,...                     % option for RSQP
    'preProcess','off',...              % option for RSQP
    'SOC',2);                           % option for RSQP

% If just 'defaults' passed in, return the default options in X
if nargin == 1 & nargout <= 1 & isequal(FUN,'defaults')
   X = defaultopt;
   return
end

medium_large = 'medium or large scale';

if nargin < 6
    error('RSQP requires at least four valid input arguments')
end
if nargin < 10, options = [];
    if nargin < 9, NONLCON=[];
        if nargin < 8, UB = [];
            if nargin < 7, LB = [];
            end, end, end, end
if isempty(NONLCON) & isempty(Aeq) & isempty(UB) & isempty(LB)
    error('RSQP is for constrained problems. Use FMINUNC for unconstrained problems.')
end

caller = 'RSQP';
XOUT = X(:);
numberOfVariables = length(XOUT);

switch optimget2(options,'Display',defaultopt,'fast')
    case {'off','none'}
        verbosity = 0;
    case 'iter'
        verbosity = 2;
    case 'final'
        verbosity = 1;
    otherwise
        verbosity = 1;
end

A = []; B = [];
% Treat with dependent linear equality constraints
if strcmp(optimget2(options,'preProcess',defaultopt,'fast'),'on')
    temp = Aeq';
    [R,jb] = rref(temp);
    temp = temp(:,jb);
    Aeq = temp';
    Beq = Beq(jb);
end

[XOUT,l,u,msg] = checkbounds(XOUT,LB,UB,numberOfVariables); 
if ~isempty(msg)
    EXITFLAG = -1;
    [FVAL,OUTPUT,LAMBDA,GRAD] = deal([]);
    X(:)=XOUT;
    if verbosity > 0
        disp(msg)
    end
    return
end

meritFunctionType = 0;
mtxmpy = optimget2(options,'HessMult',defaultopt,'fast');
if isequal(mtxmpy,'hmult')
    warnstr = sprintf('%s\n%s\n%s\n', ...
            'Potential function name clash with a Toolbox helper function:',...
            ' Use a name besides ''hmult'' for your HessMult function to',...
            '  avoid errors or unexpected results.');
    warning(warnstr)
end

diagnostics = isequal(optimget2(options,'Diagnostics',defaultopt,'fast'),'on');
funValCheck = strcmp(optimget2(options,'FunValCheck',defaultopt,'fast'),'on');
gradflag = strcmp(optimget2(options,'GradObj',defaultopt,'fast'),'on');
hessflag = strcmp(optimget2(options,'Hessian',defaultopt,'fast'),'on');
if isempty(NONLCON)
    constflag = 0;
else
    constflag = 1;
end
gradconstflag = strcmp(optimget2(options,'GradConstr',defaultopt,'fast'),'on');
line_search = 1;

if ~isempty(FUN)
    [funfcn, msg] = optimfcnchk(FUN,'RSQP',length(varargin),funValCheck,gradflag,hessflag);
else
    errmsg = sprintf('%s\n%s\n', ...
           'FUN must be a function or an inline object;', ...
           ' or, FUN may be a cell array that contains these type of objects.');
    error(errmsg)
end

if constflag
    [confcn, msg] = optimfcnchk(NONLCON,'RSQP',length(varargin),funValCheck,gradconstflag,false,1);
else
    confcn{1} = '';
end

OUTPUT.algorithm = medium_large;

if hessflag
    hessflag = 0;
    warnstr = sprintf('%s\n%s\n', ...
            'Medium-scale method is a Quasi-Newton method and does not use analytic Hessian.',...
            ' Hessian flag in options will be ignored (user-supplied Hessian will not be used).');
    warning(warnstr)
    if isequal(funfcn{1},'fungradhess')
        funfcn{1} = 'fungrad';
    elseif  isequal(funfcn{1},'fun_then_grad_then_hess')
        funfcn{1} = 'fun_then_grad';
    end 
end

if strcmp(optimget2(options,'BasisChoice',defaultopt,'fast'),'orthonormal') & ...
        strcmp(optimget2(options,'BasisChange',defaultopt,'fast'),'on')
    optimset2('BasisChange','off');
    warnstr = sprintf('%s\n%s\n', ...
            'Bases change is not used for orthonormal bases.',...
            ' BasisChange in options will be ignored.');
    warning(warnstr)
end

lenvlb = length(l);
lenvub = length(u);
CHG = 1e-7*abs(XOUT)+1e-7*ones(numberOfVariables,1);

% Ensure starting point lies within bounds
i = 1:lenvlb;
lindex = XOUT(i) < l(i);
if any(lindex),
    XOUT(lindex)=l(lindex)+1e-4; 
end
i = 1:lenvub;
uindex = XOUT(i)>u(i);
if any(uindex)
    XOUT(uindex) = u(uindex);
    CHG(uindex) = -CHG(uindex);
end
X(:) = XOUT;

% Evaluate function
GRAD = zeros(numberOfVariables,1);

switch funfcn{1}
    case 'fun'
        try
            f = feval(funfcn{3},X,varargin{:});
        catch
            errmsg = sprintf('%s\n%s\n\n%s',...
                   'RSQP cannot continue because user supplied objective function', ...
                   ' failed with the following error:', lasterr);
            error(errmsg)
        end
    case 'fungrad'
        try
            [f,GRAD(:)] = feval(funfcn{3},X,varargin{:});
        catch
            errmsg = sprintf('%s\n%s\n\n%s',...
                   'RSQP cannot continue because user supplied objective function', ...
                   ' failed with the following error:', lasterr);
            error(errmsg)
        end
    case 'fun_then_grad'
        try
            f = feval(funfcn{3},X,varargin{:});
        catch
            errmsg = sprintf('%s\n%s\n\n%s',...
                   'RSQP cannot continue because user supplied objective function', ...
                   ' failed with the following error:', lasterr);
            error(errmsg)
        end
        try
            GRAD(:) = feval(funfcn{4},X,varargin{:});
        catch
            errmsg = sprintf('%s\n%s\n\n%s',...
                   'RSQP cannot continue because user supplied objective gradient function', ...
                   ' failed with the following error:', lasterr);
            error(errmsg)
        end
    otherwise
        error('Undefined calltype in RSQP')
end

% Evaluate constraints
switch confcn{1}
    case 'fun'
        try
            [ctmp,ceqtmp] = feval(confcn{3},X,varargin{:}); 
            ceq = ceqtmp(:);
            ceqGRAD = sparse(numberOfVariables,length(ceq));
        catch
            if findstr(xlate('Too many output arguments'),lasterr)
                if isa(confcn{3},'inline')
                    errmsg = sprintf('%s%s%s\n%s\n%s\n%s', ...
                           'The inline function ',formula(confcn{3}),' representing the constraints',...
                           ' must return two outputs: the nonlinear inequality constraints and', ...
                           ' the nonlinear equality constraints.  At this time, inline objects may',...
                           ' only return one output argument: use an M-file function instead.');
                elseif isa(confcn{3},'function_handle')
                    errmsg = sprintf('%s%s%s\n%s%s', ...
                           'The constraint function ',func2str(confcn{3}),' must return two outputs:',...
                           ' the nonlinear inequality constraints and', ...
                           ' the nonlinear equality constraints.');
                else
                    errmsg = sprintf('%s%s%s\n%s%s', ...
                           'The constraint function ',confcn{3},' must return two outputs:',...
                           ' the nonlinear inequality constraints and', ...
                           ' the nonlinear equality constraints.');
                end
                error(errmsg)
            else
                errmsg = sprintf('%s\n%s\n\n%s',...
                       'RSQP cannot continue because user supplied nonlinear constraint function', ...
                       ' failed with the following error:', lasterr);
                error(errmsg)
            end
        end   
    case 'fungrad'
        try
            [ctmp,ceqtmp,cGRAD,ceqGRAD] = feval(confcn{3},X,varargin{:});
            ceq = ceqtmp(:);
        catch
            errmsg = sprintf('%s\n%s\n\n%s',...
                   'RSQP cannot continue because user supplied nonlinear constraint function', ...
                   ' failed with the following error:', lasterr);
            error(errmsg)
        end
    case 'fun_then_grad'
        try
            [ctmp,ceqtmp] = feval(confcn{3},X,varargin{:});
            ceq = ceqtmp(:);
            [cGRAD,ceqGRAD] = feval(confcn{4},X,varargin{:});
        catch
            errmsg = sprintf('%s\n%s%s\n\n%s',...
                   'RSQP cannot continue because user supplied nonlinear constraint function', ...
                   'or nonlinear constraint gradient function',...
                   ' failed with the following error:', lasterr);
            error(errmsg)
        end
    case ''
        ceq =[];
        ceqGRAD = sparse(numberOfVariables,length(ceq));
    otherwise
        error('Undefined calltype in RSQP')
end

non_eq = length(ceq);
[lin_eq,Aeqcol] = size(Aeq);
[ceqgrow, ceqgcol] = size(ceqGRAD);

eq = non_eq+lin_eq;

if ~isempty(Aeq) & Aeqcol ~= numberOfVariables
    error('Aeq has the wrong number of columns.')
end
if ceqgrow ~= numberOfVariables & ceqgcol ~= non_eq
    error('Gradient of the nonlinear equality constraints is the wrong size.')
end

% Do diagnostics on information so far
if diagnostics > 0
    msg = diagnose('rsqp',OUTPUT,gradflag,hessflag,constflag,gradconstflag,...
                   line_search,options,defaultopt,XOUT,non_eq,...
                   0,lin_eq,0,l,u,funfcn,confcn,f,GRAD,[],c,ceq,cGRAD,ceqGRAD);
end
  
% Call algorithm
[X,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN] = ...
    nlconst_rsqp(funfcn,X,l,u,Aeq,Beq,confcn,options,defaultopt, ...
                 verbosity,gradflag,gradconstflag,meritFunctionType,...
                 CHG,f,GRAD,ceq,ceqGRAD,varargin{:});

