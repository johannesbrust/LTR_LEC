function [x,FVAL,lambda_out,EXITFLAG,OUTPUT,GRADIENT,HESS] = ...
	nlconst_rsqp(funfcn,x,lb,ub,Aeq,Beq,confcn,OPTIONS,defaultopt,...
	verbosity,gradflag,gradconstflag,meritFunctionType,...
	CHG,fval,gval,nceqval,gnceqval,varargin);
%   NLCONST_RSQP Helper function to find the constrained minimum of a function 
%   of several variables. Called by RSQP.

%   Copyright 2006 Institute of Industrial Control, Zhejiang University.
%   $Revision: 1.0 $  $Date: 2006/07/08 15:35:55 $

%
%   meritFunctionType==0 for rsqp

% Initialize some parameters
FVAL = []; lambda = []; OUTPUT = []; lambdaNLP = []; HESS = [];

global OPT_STOP OPT_STEP;
global sqpopt Nlconst;
global constflag changeflag;
global numFunEvals numGradEvals numBroy numFindiff;

OPT_STEP = 1; 
OPT_STOP = 0; 
iter = 0;
XOUT = x(:);
numberOfVariables = length(XOUT);
Nlconst = 'nlconst';
bestf = Inf; 
if isempty(confcn{1})
    constflag = 0;
else
    constflag = 1;
end
stepsize = 1;
status = 0;
EXITFLAG = 1;
EXITTEXT = 'optimal solution found';

idf = Inf;
theta_SOC = 1e-5;
param1 = 20;
param2 = 10;
param4 = 100;
changeflag = 0;
IPMULT = speye(numberOfVariables);
numChBasis = 0;
numBroy = 0;
numFindiff = 0;

% Get options
tolX = optimget2(OPTIONS,'TolX',defaultopt,'fast');
tolFun = optimget2(OPTIONS,'TolFun',defaultopt,'fast');
tolCon = optimget2(OPTIONS,'TolCon',defaultopt,'fast');
DiffMinChange = optimget2(OPTIONS,'DiffMinChange',defaultopt,'fast');
DiffMaxChange = optimget2(OPTIONS,'DiffMaxChange',defaultopt,'fast');
if DiffMinChange >= DiffMaxChange
    errmsg = sprintf('%s %0.5g %s %0.5g%s\n%s\n%s',...
           'DiffMinChange options parameter is',DiffMinChange,'and DiffMaxChange is',...
           DiffMaxChange,'.','DiffMinChange must be strictly less than DiffMaxChange.');
    error(errmsg)
end
DerivativeCheck = strcmp(optimget2(OPTIONS,'DerivativeCheck',defaultopt,'fast'),'on');
maxFunEvals = optimget2(OPTIONS,'MaxFunEvals',defaultopt,'fast');
maxIter = optimget2(OPTIONS,'MaxIter',defaultopt,'fast');
if ischar(maxFunEvals)
    if isequal(lower(maxFunEvals),'100*numberofvariables')
        maxFunEvals = 100*numberOfVariables;
    else
        error('Option ''MaxFunEvals'' must be an integer value if not the default.')
    end
end
% Get special options for RSQP
BasisChange = strcmp(optimget2(OPTIONS,'BasisChange',defaultopt,'fast'),'on');
BasisChoice = optimget2(OPTIONS,'BasisChoice',defaultopt,'fast');
Correction = optimget2(OPTIONS,'Correction',defaultopt,'fast');
LSMerit = (optimget2(OPTIONS,'LSMerit',defaultopt,'fast') == 1);
SOC = optimget2(OPTIONS,'SOC',defaultopt,'fast');

sqpopt = struct('MaxSQPIter',Inf);

% Handle bounds as linear constraints
arglb = ~isinf(lb);
argub = ~isinf(ub);
if (nnz(arglb)+nnz(argub)) > 0
    A = [-speye(numberOfVariables); speye(numberOfVariables)];
    B = [-lb;ub];
    sparse(B);
else
    A = sparse(2*numberOfVariables,numberOfVariables);
    B = sparse(2*numberOfVariables,1);
end
if isempty(Aeq)
    Aeq = sparse(0,numberOfVariables); Beq = sparse(0,1);
end

x(:) = XOUT;
f = fval;
nceq = sparse(nceqval);
c = [Aeq*XOUT-Beq; nceq; A*XOUT-B];

% Get information on the number and type of constraints.
non_eq = length(nceq);
[lin_eq,Aeqcol] = size(Aeq);
[lin_ineq,Acol] = size(A); % includes upper and lower bounds
eq = non_eq+lin_eq;
ineq = lin_ineq;
ncstr = ineq+eq;

% Compute the initial constraint violation.
ga = [abs(c((1:eq)')); c((eq+1:ncstr)')];
if ~isempty(c)
    mg = max(max(ga),0);
else
    mg = 0;
end

if isempty(f)
    error('FUN must return a non-empty objective function.')
end

% Evaluate initial analytic gradients and check size.
if gradflag | gradconstflag
    if gradflag
        gf = gval;
    end
    if gradconstflag
        gnceq = gnceqval; % Don't include A and Aeq yet
        
    else
        gnceq = [];
    end
    clear nceqval gval gnceqval
end

OLDX = XOUT;
OLDgf = sparse(numberOfVariables,1);
lambdaNLP = sparse(ncstr,1);
optimError = [];

% Display header information.
header = sprintf(['\n                               max                   Directional   First-order  Hessian \n',...
        ' Iter F-count        f(x)   constraint    Step-size   derivative    optimality  Procedure ']);
formatstr = '%5.0f  %5.0f %12.6g %12.4g %12.3g %12.3g %12.3g %s  %s';
if verbosity > 1
    disp(header)
end
numFunEvals = 1;
numGradEvals = 1;
GNEW = 1e8*CHG;
how = ''; 

TIME_findiffgrad = 0;
TIME_matrixdecom = 0;
TIME_direction = 0;
TIME_qp = 0;
TIME_linesearch = 0;
TIME_broyden = 0;
TIME_findiff = 0;
TIME_bfgs = 0;

bestx = []; bsetf = []; bestgrad = []; bestlambda = [];
tic
%---------------------------------Main Loop-----------------------------
while status ~= 1
   
    %----------------GRADIENTS----------------
   
    if  ~gradconstflag | ~gradflag | DerivativeCheck
        % Finite Difference gradients (even if just checking analytical)
        TIME_findiffgrad_begin = clock;
        oldf = f;
        oldnc = nceq;
        gncFD = sparse(numberOfVariables, non_eq);
        % Try to make the finite differences equal to 1e-8.
        CHG = -1e-8./(GNEW+eps);
        % Try to make the perturbation between DiffMinChange and DiffMaxChange
        CHG = sign(CHG+eps).*min(max(abs(CHG),DiffMinChange),DiffMaxChange);    
        OPT_STEP = 1;
        for gcnt = 1:numberOfVariables
            temp = XOUT(gcnt);
            XOUT(gcnt) = temp + CHG(gcnt);
            
            % Enforce bounds while finite-differencing.
            % Need lb(gcnt) ~= ub(gcnt), and lb(gcnt) <= temp <= ub(gcnt) to enforce bounds.
            % (If the last qpsub problem was 'infeasible', the bounds could be currently violated.)
            if (lb(gcnt) ~= ub(gcnt)) & (temp >= lb(gcnt)) & (temp <= ub(gcnt)) 
                if  ((XOUT(gcnt) > ub(gcnt)) | (XOUT(gcnt) < lb(gcnt))) % outside bound ?
                    CHG(gcnt) = -CHG(gcnt);
                    XOUT(gcnt) = temp+CHG(gcnt);
                    if (XOUT(gcnt) > ub(gcnt)) | (XOUT(gcnt) < lb(gcnt)) % outside other bound ?
                        [newchg ,indsign] = max([temp-lb(gcnt), ub(gcnt)-temp]); % largest distance to bound
                        if newchg >= DiffMinChange
                            CHG(gcnt) = ((-1)^indsign)*newchg; % make sure sign is correct
                            XOUT(gcnt) = temp + CHG(gcnt);
                            warnmsg = sprintf('%s\n%s\n%s %d%s %0.5g %s\n',...
                                    'Derivative finite-differencing step was artificially reduced to be within',...
                                    ' bound constraints. This may adversely affect convergence. Increasing distance between',...
                                    '  bound constraints, in dimension',gcnt,', to be at least',abs(2*CHG(gcnt)),'may improve results.');
                            warning(warnmsg)
                        else
                            errmsg = sprintf('%s %d %s\n%s\n%s %0.5g%s\n',...
                                   'Distance between lower and upper bounds, in dimension',gcnt,', is too small to compute', ...
                                   ' finite-difference approximation of derivative. Increase distance between these', ...
                                   '  bounds to be at least',2*DiffMinChange,'.');
                            error(errmsg)
                        end
                    end
                end
            end

            x(:) = XOUT;
            if ~gradflag | DerivativeCheck
                f = feval(funfcn{3},x,varargin{:});
                gfFD(gcnt,1) = (f-oldf)/CHG(gcnt);
            end
            if ~gradconstflag | DerivativeCheck
                if constflag
                    [nc,nceq] = feval(confcn{3},x,varargin{:});
                    nceq = sparse(nceq(:));
                end
                if ~isempty(nceq)
                    gncFD(gcnt,:) = (nceq - oldnc)'/CHG(gcnt);
                end
            end
         
            OPT_STEP = 0;
            if OPT_STOP
                break;
            end
            XOUT(gcnt) = temp;
            
        end % Finish computing the finite differences of obj and constraints
        TIME_findiffgrad = TIME_findiffgrad+etime(clock,TIME_findiffgrad_begin);
        % Adjust to new order after finite differences
        if changeflag
            XOUT = IPMULT'*XOUT;
            [reorder,j] = find(IPMULT);
            lb = lb(reorder);
            ub = ub(reorder);
        end
        
        % Gradient check
        if DerivativeCheck == 1 & (gradflag | gradconstflag) % Analytic exists                               
            disp('Function derivative')
            if gradflag                
                if isa(funfcn{4},'inline')
                    graderr(gfFD, gf, formula(funfcn{4}));
                else
                    graderr(gfFD, gf, funfcn{4});
                end
            end
         
            if gradconstflag
                disp('Constraint derivative')
                if isa(confcn{4},'inline')
                    graderr(gncFD, gnceq, formula(confcn{4})); 
                else
                    graderr(gncFD, gnceq, confcn{4});
                end
            end         
            DerivativeCheck = 0;
        end % DerivativeCheck == 1 & (gradflag | gradconstflag)
        
        if ~gradflag
            gf = gfFD;
        end
        if ~gradconstflag
            gnceq = gncFD;
        end
        
        numFunEvals = numFunEvals + numberOfVariables;
        f = oldf;
        nceq = oldnc;
        clear oldf oldnc
    end
    
    if changeflag
        gf = IPMULT'*gf;
        % Now add in Aeq and A
        if ~isempty(gnceq)
            gc = IPMULT'*[Aeq',gnceq,A'];
        elseif ~isempty(Aeq) | ~isempty(A)
            gc = IPMULT'*[Aeq',A'];
        else
            gc = sparse(numberOfVariables,0);
        end
    else
        % Now add in Aeq and A
        if ~isempty(gnceq)
            gc = [Aeq',gnceq,A'];
        elseif ~isempty(Aeq) | ~isempty(A)
            gc = [Aeq',A'];
        else
            gc = sparse(numberOfVariables,0);
        end
    end
    
    OPT_STEP = 2;
    
    if iter > 0
       
        % Compute the first order KKT conditions.       
        normgradLag = norm(gf + gc*lambdaNLP,inf);
        temp = c(eq+1:ncstr);
        temp(isinf(temp)) = 0;
        normcomp = norm(lambdaNLP(eq+1:ncstr).*temp,inf);
        clear temp;
        if isfinite(normgradLag) & isfinite(normcomp)
            optimError = max(normgradLag, normcomp);
        else
            optimError = inf;
        end
        feasError = mg;
        optimScal = 1; feasScal = 1;
        
        if changeflag
            gfSD = (IPMULT*gf)'*SD;
        else
            gfSD = gf'*SD;
        end
            
        % Print iteration information     
        if verbosity > 1 
            CurrOutput = sprintf(formatstr,iter,numFunEvals,f,full(mg),...
                                 stepsize,gfSD,optimError,how,howqp); 
            disp(CurrOutput)
        end
        
        %-------------TEST CONVERGENCE---------------

        if optimError < tolFun*optimScal & feasError < tolCon*feasScal
            if verbosity > 0
                disp('Optimization terminated successfully:')
                disp(' First-order optimality measure less than OPTIONS.TolFun and')
                disp('  maximum constraint violation is less than OPTIONS.TolCon')
            end
            EXITFLAG = 1;
            status = 1;
            active_const = find(lambda ~= 0);
            if active_const 
                if verbosity > 0
                    if size(active_const,1) < 100
                        disp('Active Constraints:')
                        disp(active_const') 
                    else
                        fprintf('Number of Active Constraints: %d\n',size(active_const,1))
                    end
                end      
            else % active_const == 0
                if verbosity > 0
                    disp(' No Active Constraints')
                end
            end % active_const
        elseif (max(abs(SD)) < 2*tolX | abs(gfSD) < 2*tolFun) & ...
               (mg < tolCon | (strncmp(howqp,'i',1) & mg > 0))
            % The algorithm can make no more progress.  If feasible, compute 
            % the new up-to-date Lagrange multipliers (with new gradients) 
            % and recompute the KKT error.  Then output appropriate termination
            % message.
            if ~strncmp(howqp, 'i', 1)
                if optimError < tolFun*optimScal
                    if verbosity > 0
                        disp('Optimization terminated successfully:')
                        disp(' First-order optimality measure less than OPTIONS.TolFun and')
                        disp('  maximum constraint violation is less than OPTIONS.TolCon')
                    end
                    EXITFLAG = 1;
                elseif max(abs(SD)) < 2*tolX
                    if verbosity > 0
                        disp('Optimization terminated successfully:')
                        disp(' Search direction less than 2*OPTIONS.TolX and')
                        disp('  maximum constraint violation is less than OPTIONS.TolCon')
                    end
                    EXITFLAG = 4;
                else 
                    if verbosity > 0 
                        disp('Optimization terminated successfully:')
                        disp(' Magnitude of directional derivative in search direction ')
                        disp('  less than 2*OPTIONS.TolFun and maximum constraint violation ')
                        disp('   is less than OPTIONS.TolCon') 
                    end
                    EXITFLAG = 5;
                end  
                
                active_const = find(lambda ~= 0);
                if active_const 
                    if verbosity > 0
                        if size(active_const,1) < 100
                            disp('Active Constraints:') 
                            disp(active_const')
                        else
                            fprintf('Number of Active Constraints: %d\n',size(active_const,1))
                        end
                    end
                else
                    if verbosity > 0
                        disp('No Active Constraints')
                    end   
                end

            end % ~strncmp(howqp, 'i', 1)
        
            if (strncmp(howqp, 'i',1) & mg > 0)
                if verbosity > 0
                    disp('Optimization terminated: No feasible solution found.')
                end 
                if max(abs(SD)) < 2*tolX
                    if verbosity > 0
                        disp(' Search direction less than 2*OPTIONS.TolX but constraints are not satisfied.')
                    end
                else
                    if verbosity > 0
                        disp(' Magnitude of directional derivative in search direction ')
                        disp('  less than 2*OPTIONS.TolFun but constraints are not satisfied.')
                    end
                end
                EXITFLAG = -2;  
                EXITTEXT = 'no feasible solution';
            end
            status = 1;
        else % continue
            if numFunEvals > maxFunEvals  | OPT_STOP
                XOUT = OLDX;
                f = OLDF;
                if ~OPT_STOP
                    if verbosity > 0
                        disp('Maximum number of function evaluations exceeded;')
                        disp('increase OPTIONS.MaxFunEvals')
                    end
                end
                EXITFLAG = 0;
                EXITTEXT = 'maxFunEvals exceeded';
                status = 1;
            end
            if iter > maxIter
                XOUT = OLDX;
                f = OLDF;
                if verbosity > 0
                    disp('Maximum number of iterations exceeded;')
                    disp('increase OPTIONS.MaxIter')
                end
                EXITFLAG = 0;
                EXITTEXT = 'maxIter exceeded';
                status = 1;
            end
        end
    end
    if status ~= 1
        how= ''; 
        iter = iter + 1; 
        
        if numGradEvals > 1 % Check for first call
            OLDY = Y;
            OLDZ = Z;
            [Y,Z,gcY,P,NewBasis,idf,curTIME_matrixdecom] = ...
                SPACE_DECOM(gc(:,1:eq),stepsize,numChBasis,idf,BasisChange,BasisChoice);
            if NewBasis
                changeflag = 1;
                numChBasis = numChBasis+1;
                if verbosity > 0
                    fprintf('Basis changed %d time(s).\n',numChBasis)
                end
                IPMULT = IPMULT*P';
                XOUT = P*XOUT;
                gf = P*gf;
                gc = P*gc;
                OLDX = P*OLDX;
                OLDgf = P*OLDgf;
                OLDgc = P*OLDgc;
                [reorder,j] = find(P');
                lb = lb(reorder);
                ub = ub(reorder);
                lambda_bnd = P*lambda_bnd;
                temp = sparse(2*numberOfVariables,1);
                temp(ACTIND_lu) = 1;
                temp = [P*temp(1:numberOfVariables); P*temp(numberOfVariables+1:end)];
                ACTIND_lu = find(temp);
                ACTIND_lu = ACTIND_lu';
                clear temp
            end
            
        else % First call
            lambda_lu = sparse(2*numberOfVariables,1);
            lambda_bnd = sparse(numberOfVariables,1);
            ACTIND_eq = 1:eq;
            ACTIND_lu = [];
            
            [Y,Z,gcY,P,NewBasis,idf,curTIME_matrixdecom] = ...
                SPACE_DECOM(gc(:,1:eq),stepsize,numChBasis,idf,BasisChange,BasisChoice);
            if NewBasis
                changeflag = 1;
                numChBasis = numChBasis+1; 
                if verbosity > 0
                    fprintf('Basis changed %d time(s).\n',numChBasis)
                end
                IPMULT = IPMULT*P'; 
                XOUT = P*XOUT;
                gf = P*gf;
                gc = P*gc;
                OLDX = P*OLDX;
                OLDgf = P*OLDgf;
                [reorder,j] = find(P');
                lb = lb(reorder);
                ub = ub(reorder);
            end
            
            OLDLAMBDA = (eps+gf'*gf)*ones(ncstr,1)./(sum(gc.*gc)'+eps);
            
        end % if numGradEvals > 1
        
        %-------------SEARCH DIRECTION---------------
        TIME_direction_begin = clock;       
        
        if strcmp(BasisChoice,'coordinate')
            temp = Y'*(gf+lambda_bnd);
            lambda_eq = -gcY'\temp;
        else
            temp = Y'*(gf+lambda_bnd);
            lambda_eq = -(Y'*gc(:,1:eq))\temp;
        end
        clear temp;
        
        if numGradEvals > 1
            % Update Broy_A and RHess
            if NewBasis
                Broy_A = [sparse(numberOfVariables-eq,eq),speye(numberOfVariables-eq)];  
                RHess = speye(numberOfVariables-eq);
                curTIME_broyden_rh = 0;
                curTIME_findiff_rh = 0;
                curTIME_bfgs = 0;
            else                
                if strcmp(BasisChoice,'coordinate')
                    broy_y = OLDZ'*(gf+gc(:,1:eq)*lambda_eq-OLDgf);
                else
                    broy_y = OLDZ'*(gf-OLDgf+(gc(:,1:eq)-OLDgc(:,1:eq))*lambda_eq);
                end
                broy_s = XOUT-OLDX;
                Broy_A = Broy_A+(broy_y-Broy_A*broy_s)*broy_s'/(broy_s'*broy_s);
                
                [RHess,curTIME_broyden_rh,curTIME_findiff_rh,curTIME_bfgs,how] = ...
                    UPDATE_RHESS(funfcn,XOUT,OLDX,Aeq,confcn,BasisChoice,Correction,OLDF,OLDgf,OLDgc(:,lin_eq+1:eq),Broy_A, ...
                	RHess,broy_y,OLDY,OLDZ,py,pz,lambda_eq,IPMULT,findiff,regionflag,stepsize,param5,param6,varargin{:});
            end
        else
            % Initialize Broy_A and RHess
            Broy_A = [sparse(numberOfVariables-eq,eq),speye(numberOfVariables-eq)];
            RHess = speye(numberOfVariables-eq);
            curTIME_broyden_rh = 0;
            curTIME_findiff_rh = 0;
            curTIME_bfgs = 0;
        end
        
        % step in range space
        if strcmp(BasisChoice,'coordinate')
            py = -gcY\c(1:eq);
        else
            py = -(gc(:,1:eq)'*Y)\c(1:eq);
        end
        Ypy = Y*py;
        
        param3 = 0.1*(numberOfVariables-eq)^0.25/iter^1.1;
        param5 = 0.1*(numberOfVariables-eq)^0.25/iter^1.1;
        param6 = 0.01*(numberOfVariables-eq)^0.25/iter^1.1;
        % step in null space
        [pz,w,pdamp,findiff,regionflag,lambda_lu,lambda_bnd,ACTIND_lu,howqp,curTIME_broyden_pz,curTIME_findiff_pz,curTIME_qp] = ...
            GET_PZ(funfcn,XOUT,lb,ub,A,B,Aeq,confcn,OPTIONS,BasisChoice,Correction,gf,c(1:eq),gc(:,lin_eq+1:eq),Broy_A,RHess, ...
            Ypy,Z,py,lambda_eq,lambda_bnd,ACTIND_lu,IPMULT,param1,param2,param3,param4,gc(:,1:eq),varargin{:});
        if changeflag
            lambda_lu = [IPMULT*lambda_lu(1:numberOfVariables); IPMULT*lambda_lu(numberOfVariables+1:end)];
        end
        Zpz = Z*pz;
        
        % search direction
        SD = Ypy+Zpz;
        TIME_direction = TIME_direction+etime(clock,TIME_direction_begin);
        
        lambda = [lambda_eq; lambda_lu];
        if changeflag
            temp = sparse(2*numberOfVariables,1);
            temp(ACTIND_lu) = 1;
            temp = [IPMULT*temp(1:numberOfVariables); IPMULT*temp(numberOfVariables+1:end)];
            temp = find(temp);
            temp = temp';
        else
            temp = ACTIND_lu;
        end
        ACTIND = [ACTIND_eq,temp+eq];
        clear temp
        lambdaNLP(:,1) = 0;
        lambdaNLP(ACTIND) = lambda(ACTIND);
        % penalty factor
        if  norm(OLDLAMBDA,inf) < 1.001*norm(lambda,inf)
            OLDLAMBDA = 1.001*[abs(lambda(1:eq)); lambda(eq+1:end)];
        else    
            OLDLAMBDA = 0.25*(3*OLDLAMBDA+[abs(lambda(1:eq)); lambda(eq+1:end)]);
        end
        
        ga = [abs(c(1:eq)); c(eq+1:ncstr)];
        if ~isempty(c)
            mg = max(max(ga),0);
        else
            mg = 0;
        end
        
        if strncmp(howqp,'ok',2); 
            howqp = ''; 
        end
        if ~isempty(how) & ~isempty(howqp) 
            how = [how,'; '];
        end
        
        OLDX = XOUT;
        OLDF = f;
        OLDgf = gf;
        OLDgc = gc;
        XN = sparse(numberOfVariables,1);
        
        %---------------LINESEARCH--------------------
        TIME_linesearch_begin = clock;

        temp = ga;
        temp(isinf(ga) | (ga < 0)) = 0;
        MATL = f+sum(OLDLAMBDA.*temp);
        
        if LSMerit
        	DGMATL = gf'*SD-sum(OLDLAMBDA.*temp);
        else        
            infeas = strncmp(howqp,'i',1);
            if mg > 0
                MATL2 = mg;
            elseif f >= 0 
                MATL2 = -1/(f+1);
            else 
                MATL2 = 0;
            end
            if ~infeas & f < 0
                MATL2 = MATL2 + f - 1;
            end
        end
        
        if changeflag
            MATX = IPMULT*XOUT;
            SD = IPMULT*SD;
        else
            MATX = XOUT;
        end
        
        if mg < eps & f < bestf
            bestf = f;
            if changeflag 
                bestx = IPMULT*XOUT;
                bestgrad = IPMULT*gf;
            else
                bestx = XOUT;
                bestgrad = gf;
            end
            bestlambda = lambda;
        end
        
        x = MATX+SD;
        f = feval(funfcn{3},x,varargin{:});
        numFunEvals = numFunEvals + 1;
        if isnan(f)
            EXITFLAG = -4;
            EXITTEXT = 'NaN occurred';
            status = 1;
            break;
        end
        
        if constflag
            [nc,nceq] = feval(confcn{3},x,varargin{:});
        else
            nceq = [];
        end
        c = [Aeq*x-Beq; nceq(:); A*x-B];
        ga = [abs(c(1:eq)); c(eq+1:ncstr)];
        if ~isempty(c)
            mg = max(max(ga),0);
        else
            mg = 0;
        end
        
        temp = ga;
        temp(isinf(ga) | (ga < 0)) = 0;
        MERIT = f+sum(OLDLAMBDA.*temp);
        if LSMerit
            MeritCondition = (MERIT > MATL+0.1*DGMATL);
        else
            if mg > 0
               MERIT2 = mg;
            elseif f >=0
               MERIT2 = -1/(f+1);
            else
               MERIT2 = 0;
            end
            if ~infeas & f < 0
               MERIT2 = MERIT2 + f - 1;
            end
            MeritCondition = ((MERIT2 > MATL2) & (MERIT > MATL));
        end
        
        % Treat with Maratos effect
        SOC_flag = 0;
        if MeritCondition & (mg < theta_SOC)
            if LSMerit
                param = DGMATL;
                infeas = [];
            else
                param = MATL2;
            end
            if SOC == 1
                if strcmp(BasisChoice,'coordinate')
                    gceqY = gcY;
                else
                    gceqY = gc(:,1:eq)'*Y;
                end
                [x,c,mg,lambda_lu_c,ACTIND_lu,SOC_flag] = ...
                    TRY_SOC(funfcn,XOUT,lb,ub,A,B,Aeq,Beq,confcn,OPTIONS,gf,ga,RHess,Ypy,Zpz,Y,Z,pz, ...
                    MATL,param,w,pdamp,eq,OLDLAMBDA,ACTIND_lu,IPMULT,infeas,gceqY,c(1:eq),Broy_A,SOC,varargin{:});
            else
                [x,c,mg,lambda_lu_c,ACTIND_lu,SOC_flag] = ...
                    TRY_SOC(funfcn,XOUT,lb,ub,A,B,Aeq,Beq,confcn,OPTIONS,gf,ga,RHess,Ypy,Zpz,[],Z,pz, ...
                    MATL,param,w,pdamp,eq,OLDLAMBDA,ACTIND_lu,IPMULT,infeas,[],[],[],SOC,varargin{:});
            end
            if SOC_flag > 0
                lambda_bnd = lambda_lu_c(numberOfVariables+1:2*numberOfVariables,1)-lambda_lu_c(1:numberOfVariables,1);
                if changeflag
                    lambda_lu = [IPMULT*lambda_lu_c(1:numberOfVariables); IPMULT*lambda_lu_c(numberOfVariables+1:end)];
                else
                    lambda_lu = lambda_lu_c;
                end
            end
        end
        
        stepsize = 1;
        while MeritCondition & (numFunEvals < maxFunEvals) & ~OPT_STOP & (SOC_flag <= 0)
            
        	stepsize = stepsize/2;
            if stepsize < 1e-4
                stepsize = -stepsize; 
            end
            
            x = MATX+stepsize*SD;
            f = feval(funfcn{3},x,varargin{:});
            numFunEvals = numFunEvals + 1;
            if isnan(f)
                EXITFLAG = -4;
                OUTPUT.ExitText = 'NaN occurred';
                status = 1;
                break;
            end
            
            if constflag
                [nc,nceq] = feval(confcn{3},x,varargin{:});
            else
                nceq = [];
            end
            c = [Aeq*x-Beq; nceq(:); A*x-B];
            ga = [abs(c(1:eq)); c(eq+1:ncstr)];
            if ~isempty(c)
                mg = max(max(ga),0);
            else
                mg = 0;
            end
            
            if OPT_STOP
                break;
            end
            
            temp = ga;
            temp(isinf(ga) | (ga < 0)) = 0;
            MERIT = f+sum(OLDLAMBDA.*temp);
            if LSMerit
                MeritCondition = (MERIT > MATL+0.1*stepsize*DGMATL);
            else
                if mg > 0
                   MERIT2 = mg;
                elseif f >=0 
                   MERIT2 = -1/(f+1);
                else 
                   MERIT2 = 0;
                end
                if ~infeas & f < 0
                   MERIT2 = MERIT2 + f - 1;
                end
                MeritCondition = ((MERIT2 > MATL2) & (MERIT > MATL));
            end
            
        end  % line search loop
        clear temp;
        
        TIME_linesearch = TIME_linesearch+etime(clock,TIME_linesearch_begin);
        TIME_matrixdecom = TIME_matrixdecom+curTIME_matrixdecom;
        TIME_broyden = TIME_broyden+curTIME_broyden_rh+curTIME_broyden_pz;
        TIME_findiff = TIME_findiff+curTIME_findiff_rh+curTIME_findiff_pz;
        TIME_qp = TIME_qp+curTIME_qp;
        TIME_bfgs = TIME_bfgs+curTIME_bfgs;
        
        if EXITFLAG == -4
            break
        end
        
        if ~gradconstflag | ~gradflag
            % Need to be in original order until finishing calculation of finite differences
            XOUT = x;
            if changeflag
                [reorder,j] = find(IPMULT');
                lb = lb(reorder);
                ub = ub(reorder);
            end
        elseif changeflag
            XOUT = IPMULT'*x;
        else
            XOUT = x;
        end       
        
        % Evaluate function gradients
        switch funfcn{1}
            case 'fun'
                ;   % do nothing...will use finite difference
            case 'fungrad'
                [f,gf] = feval(funfcn{3},x,varargin{:});
                numGradEvals = numGradEvals+1;
            case 'fun_then_grad'
                gf = feval(funfcn{4},x,varargin{:});
                numGradEvals = numGradEvals+1;
            otherwise
                error('Undefined calltype in RSQP')
        end
        numFunEvals = numFunEvals+1;
   
        % Evaluate constraint gradients
        switch confcn{1}
            case 'fun'
                gnceq = [];
            case 'fungrad' 
                [nc,nceq,gncineq,gnceq] = feval(confcn{3},x,varargin{:});
                numGradEvals = numGradEvals+1;
            case 'fun_then_grad'
                [gncineq,gnceq] = feval(confcn{4},x,varargin{:});
                numGradEvals = numGradEvals+1;
            case ''
                gnceq = sparse(numberOfVariables,0);
            otherwise
                error('Undefined calltype in RSQP')
        end
        
        if toc >= 10800
            EXITFLAG = -3;
            EXITTEXT = 'CPU time limit was exceeded';
            status = 1;
        end
    end % if status ~= 1
end % while status ~= 1

% Update 
numConstrEvals = numGradEvals;

% If a better unconstrained solution was found earlier, use it:
if (f > bestf) | isnan(f) 
	x = bestx;
	f = bestf;
	GRADIENT = bestgrad;
	lambda = bestlambda;
elseif changeflag
    x(:) = IPMULT*XOUT;
	GRADIENT = IPMULT*gf;
else
    x(:) =XOUT;    
    GRADIENT = gf; 
end

FVAL = f;

if (OPT_STOP)
    EXITFLAG = -1;
    EXITTEXT = 'Optimization terminated prematurely by user';
    if verbosity > 0
        disp('Optimization terminated prematurely by user')
    end
end

OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.stepsize = stepsize;
OUTPUT.algorithm = 'medium-large scale: RSQP, Quasi-Newton, line-search';
OUTPUT.firstorderopt = optimError;
OUTPUT.cgiterations = [];

lambda_out.lower = sparse(numberOfVariables,1);
lambda_out.upper = sparse(numberOfVariables,1);
lambda_out.eqlin = lambdaNLP(1:lin_eq);
ii = lin_eq ;
lambda_out.eqnonlin = lambdaNLP(ii+1: ii+non_eq);
ii = ii+non_eq;
lambda_out.lower = lambdaNLP(ii+1:ii+numberOfVariables);
ii = ii + numberOfVariables;
lambda_out.upper = lambdaNLP(ii+1:end);

OUTPUT.FVAL = FVAL;
OUTPUT.x = x;
OUTPUT.ExitFlag = EXITFLAG;
OUTPUT.ExitText = EXITTEXT;
OUTPUT.numChBasis = numChBasis;
OUTPUT.numBroy = numBroy;
OUTPUT.numFindiff = numFindiff;
OUTPUT.UsedTime.findiffgrad = TIME_findiffgrad;
OUTPUT.UsedTime.spacedecom = TIME_matrixdecom;
OUTPUT.UsedTime.direction = TIME_direction;
OUTPUT.UsedTime.qp = TIME_qp;
OUTPUT.UsedTime.linesearch = TIME_linesearch;
OUTPUT.UsedTime.broyden = TIME_broyden;
OUTPUT.UsedTime.findiff = TIME_findiff;
OUTPUT.UsedTime.bfgs = TIME_bfgs;

% NLCONST_RSQP finished

%--------------------------------------------------------------------------

function [Y,Z,gcY,P,NewBasis,idf,Time_matrixdecom] = ...
    SPACE_DECOM(gceq,stepsize,numChBasis,idf,BasisChange,BasisChoice)
%Implement space decompositon and basis change if necessary.

Time_matrixdecom_begin = clock;

gceq = gceq';
[eq,numberOfVariables] = size(gceq);

if sprank(gceq) < eq
    error('Dependent constraints cannot be decomposed')
end

P = [];
NewBasis = 0;

if strcmp(BasisChoice,'orthonormal')
    gcY = [];
else
    gcY = gceq(:,1:eq);
    gcX = gceq(:,eq+1:numberOfVariables);
    gcYX = gcY\gcX;
    if BasisChange
        idf_now = max(max(abs(gcYX)));
        if  sprank(gcY) < eq |...
            nnz(isinf(gcYX)) | nnz(isnan(gcYX)) | ...
            ((stepsize < 1e-3/2^numChBasis) & (abs(idf_now) > 2*numChBasis*abs(idf))) 
            % find nonsingular leading square matrix through the interface
            % of optim toolbox for matlab
            %result = proxy('findp',1,gceq');
            % 06/21/18, J.B. Modification of proxy call
            %result = proxy('findp',1,gceq');
            %P = result{1};
            P = findp(gceq');
            
            gceq = gceq*P';
            gcY = gceq(:,1:eq);
            if sprank(gcY) < eq
                error('Dependent constraints cannot be decomposed')
            end
            gcX = gceq(:,eq+1:numberOfVariables);
            gcYX = gcY\gcX;
            idf = max(max(gcYX));
            NewBasis = 1;
        else
            idf = idf_now;
        end
    end       
end

switch BasisChoice
    case 'orthonormal'
        [Q R] = qr(gceq');
        Y = Q(1:numberOfVariables,1:eq);
        Z = Q(1:numberOfVariables,eq+1:numberOfVariables); 
        sparse(Y);
        sparse(Z);
    case 'orthogonal'
        Z = [-gcYX; speye(numberOfVariables-eq)];
        Y = [speye(eq); gcYX'];
        sparse(Y);
        sparse(Z);
    case 'coordinate'
        Z = [-gcYX; speye(numberOfVariables-eq)];
        Y = [speye(eq); sparse(numberOfVariables-eq,eq)];
        sparse(Z);
        sparse(Y);
end

Time_matrixdecom = etime(clock,Time_matrixdecom_begin);

%--------------------------------------------------------------------------

function [pz,w,pdamp,findiff,regionflag,lambda_lu,lambda_bnd,ACTIND_lu,howqp,Time_broyden,Time_findiff,Time_qp] = ...
	GET_PZ(funfcn,XOUT,lb,ub,A,B,Aeq,confcn,OPTIONS,BasisChoice,Correction,gf,ceqval,gnceq,Broy_A,RHess, ...
	Ypy,Z,py,lambda_eq,lambda_bnd,ACTIND_lu,IPMULT,param1,param2,param3,param4,gceq,varargin)
%Calculate step in null space, pz.

global Nlconst sqpopt changeflag;
global numFunEvals numGradEvals numBroy numFindiff;

numberOfVariables = length(XOUT);
w = 0;
pdamp = 0;
findiff = 0;
Time_broyden = 0;
Time_findiff = 0;

if Correction == 2  
    % compute the cross term w
    Time_broyden_begin = clock;
    w = Broy_A*Ypy;
    if norm(w) > param1*(norm(py))^0.5
        w = w*param1*(norm(py))^0.5/norm(w);
    end
    Time_broyden = etime(clock,Time_broyden_begin);

    % compute damping parameter
    pdamp=gf'*Z*(RHess\w);
    if pdamp >= 0 
        pdamp = 1;
    else
        temp = Z'*gf;
        pdamp = min(0.1*gf'*Z*(RHess\temp)/pdamp,1);    
    end
end

% compute pz
Time_qp_begin = clock;
qp_f = Z'*gf+pdamp*w; 
qp_A = [-Z;Z];
qp_B = [-lb+XOUT+Ypy; ub-XOUT-Ypy];
qp_XN = sparse(size(Z,2),1);
[pz,lambda_lu,exitflagqp,outputqp,howqp,ACTIND_lu_new] = ...
    qpsub2(RHess,qp_f,qp_A,qp_B,[],[],qp_XN,0,-1,Nlconst, ...
           size(qp_A,1),size(Z,2),OPTIONS,sqpopt,ACTIND_lu);
Time_qp = etime(clock,Time_qp_begin);

if changeflag
    x = IPMULT*XOUT;
else
    x = XOUT;
end

% Determine region of the update criterion for BFGS
opterr = norm(Z'*(gf+lambda_bnd))+norm(ceqval);
regionflag = (norm(py) <= param2*norm(pz)/opterr^0.5);
if  regionflag & (norm(py) > param3^2*norm(pz)) & ...
    (strcmp(funfcn{1},'fungrad') |  strcmp(funfcn{1},'fun_then_grad')) & ...
    (strcmp(confcn{1},'fungrad') |  strcmp(confcn{1},'fun_then_grad'))  
    findiff = 1;
    numFindiff = numFindiff+1;
else
    numBroy = numBroy+1;
end

if (Correction == 2) & findiff 
    % re-compute w via finite differences
    Time_broyden = 0;
    Time_findiff_begin = clock;
    if changeflag
        x = IPMULT*(XOUT+Ypy);
    else
        x = XOUT+Ypy;
    end

    switch funfcn{1}
        case 'fungrad'
            [f_tmp,gf_tmp] = feval(funfcn{3},x,varargin{:});
            if changeflag, gf_tmp = IPMULT'*gf_tmp(:); end
            numGradEvals = numGradEvals+1;
        case 'fun_then_grad'
            gf_tmp = feval(funfcn{4},x,varargin{:});
            if changeflag, gf_tmp = IPMULT'*gf_tmp(:); end
            numGradEvals = numGradEvals+1;
    end

    switch confcn{1}
        case 'fungrad'
            [ncineq_tmp,nceq_tmp,gncineq_tmp,gnceq_tmp] = feval(confcn{3},x,varargin{:});
            if changeflag
                gnceq_tmp = IPMULT'*gnceq_tmp;
            end
            numGradEvals = numGradEvals+1;
        case 'fun_then_grad'
            [gncineq_tmp,gnceq_tmp] = feval(confcn{4},x,varargin{:});
            if changeflag
                gnceq_tmp = IPMULT'*gnceq_tmp;              
            end
            numGradEvals = numGradEvals+1;
        case ''
            gnceq_tmp = sparse(numberOfVariables,0);
    end

    if strcmp(BasisChoice,'coordinate')
        if changeflag, Aeq = Aeq*IPMULT; end
        w = Z'*(gf_tmp+[Aeq',gnceq_tmp]*lambda_eq-gf);
    else
        lin_eq = size(Aeq,1);
        w = Z'*(gf_tmp-gf+(gnceq_tmp-gnceq)*lambda_eq(lin_eq+1:end));
    end
    if norm(w) > param4*norm(py)^0.5
        w = w*param4*norm(py)^0.5/norm(w);
    end
    Time_findiff = etime(clock,Time_findiff_begin);

    % re-compute damping parameter
    pdamp=gf'*Z*(RHess\w);
    if pdamp >= 0 
        pdamp = 1;
    else
        temp = Z'*gf;
        pdamp = min(0.1*gf'*Z*(RHess\temp)/pdamp,1);    
    end

    % re-compute pz
    Time_qp_begin = clock;
    qp_f = Z'*gf+pdamp*w; 
    qp_A = [-Z;Z];
    qp_B = [-lb+XOUT+Ypy; ub-XOUT-Ypy];
    qp_XN = sparse(size(Z,2),1);
    [pz,lambda_lu,exitflagqp,outputqp,howqp,ACTIND_lu_new] = ...
        qpsub2(RHess,qp_f,qp_A,qp_B,[],[],qp_XN,0,-1,Nlconst, ...
               size(qp_A,1),size(Z,2),OPTIONS,sqpopt,ACTIND_lu);            
    Time_qp = Time_qp+etime(clock,Time_qp_begin);
end

lambda_bnd = lambda_lu(numberOfVariables+1:end)-lambda_lu(1:numberOfVariables);
ACTIND_lu = ACTIND_lu_new;

%--------------------------------------------------------------------------

function [RHess,Time_broyden,Time_findiff,Time_bfgs,how] = ...
	UPDATE_RHESS(funfcn,XOUT,OLDX,Aeq,confcn,BasisChoice,Correction,OLDF,OLDgf,OLDgnceq,Broy_A,RHess, ...
	broy_y,OLDY,OLDZ,py,pz,lambda_eq,IPMULT,findiff,regionflag,stepsize,param5,param6,varargin)
%Updates reduced hessian with damped BFGS methods.

global changeflag numFunEvals numGradEvals numBroy numFindiff;

Time_broyden = 0;
Time_findiff = 0;

if ~Correction
    updateB_w = 0;
% Calculate the cross term updateB_w for BFGS update
elseif ~findiff
    Time_broyden_begin = clock;
    updateB_w = stepsize*Broy_A*(OLDY*py);
    if norm(updateB_w) > abs(stepsize)*norm(py)/param5
        updateB_w = updateB_w*abs(stepsize)*norm(py)/(param5*norm(updateB_w));
    end
    Time_broyden = etime(clock,Time_broyden_begin);
else
    Time_findiff_begin = clock;
    if changeflag
        x = IPMULT*(OLDX+OLDY*py);
    else
        x = OLDX+OLDY*py;
    end

    switch funfcn{1}
        case 'fungrad'
            [f_tmp,gf_tmp] = feval(funfcn{3},x,varargin{:});
            if changeflag, gf_tmp = IPMULT'*gf_tmp(:); end
            numGradEvals = numGradEvals+1;
        case 'fun_then_grad'
            gf_tmp = feval(funfcn{4},x,varargin{:});
            if changeflag, gf_tmp = IPMULT'*gf_tmp(:); end
            numGradEvals = numGradEvals+1;
    end

    switch confcn{1}
        case 'fungrad'
            [ncineq_tmp,nceq_tmp,gncineq_tmp,gnceq_tmp] = feval(confcn{3},x,varargin{:});
            if changeflag, gnceq_tmp = IPMULT'*gnceq_tmp; end           
            numGradEvals = numGradEvals+1;
        case 'fun_then_grad'
            [gncineq_tmp,gnceq_tmp] = feval(confcn{4},x,varargin{:});
            if changeflag, gnceq_tmp = IPMULT'*gnceq_tmp; end            
            numGradEvals = numGradEvals+1;
        case ''
            gnceq_tmp = sparse(numberOfVariables,0);
    end

    if strcmp(BasisChoice,'coordinate')
        if changeflag, Aeq = Aeq*IPMULT; end
        updateB_w = stepsize*OLDZ'*(gf_tmp+[Aeq',gnceq_tmp]*lambda_eq-OLDgf); 
    else
        lin_eq = size(Aeq,1);
        updateB_w = stepsize*OLDZ'*(gf_tmp-OLDgf+(gnceq_tmp-OLDgnceq)*lambda_eq(lin_eq+1:end)); % 0605
    end
    if norm(updateB_w) > abs(stepsize)*norm(py)/param6
        updateB_w = updateB_w*abs(stepsize)*norm(py)/(param6*norm(updateB_w));
    end
    Time_findiff = etime(clock,Time_findiff_begin);
end

% Update reduced hessian
Time_bfgs_begin = clock;
B_s = stepsize*pz;
B_y = broy_y-updateB_w;
if B_s'*B_y >= 0.2*B_s'*RHess*B_s
    theta = 1;
else
    theta = 0.8*B_s'*RHess*B_s/(B_s'*RHess*B_s-B_s'*B_y);
end
B_y = theta*B_y+(1-theta)*RHess*B_s;
if (B_s'*B_y) > eps & regionflag
    RHess = RHess-((RHess*B_s)*(B_s'*RHess))/(B_s'*RHess*B_s)+B_y*B_y'/(B_y'*B_s);
    how = ' Hessian updated';
else
    how = ' Hessian not updated';
end
Time_bfgs = etime(clock,Time_bfgs_begin);

%--------------------------------------------------------------------------

function [x,c,mg,lambda_lu_c,ACTIND_lu,SOC_flag] = ...
    TRY_SOC(funfcn,XOUT,lb,ub,A,B,Aeq,Beq,confcn,OPTIONS,gf,ga,RHess,Ypy,Zpz,Y,Z,pz, ...
    MATL,param,w,pdamp,eq,OLDLAMBDA,ACTIND_lu,IPMULT,infeas,gceqY,ceq_tmp,Broy_A,SOC,varargin)
%Calculate second order correction step to overcome Maratos effect.

global  Nlconst sqpopt changeflag constflag;

if SOC == 1
    % Compute the correction step in range space, py_c
    py_c = -gceqY\ceq_tmp;
    % Compute the correction step in null space, pz_c
    qp_f = Z'*gf+pdamp*w+RHess*pz+Broy_A*(Y*py_c);
    qp_A = [-Z;Z];
    qp_B = [-lb+XOUT+Ypy+Zpz+Y*py_c;ub-XOUT-Ypy-Zpz-Y*py_c];
    qp_XN = sparse(size(Z,2),1);
    [pz_c,lambda_lu_c,exitflagqp,outputqp,howqp,ACTIND_lu_c] = ...
        qpsub2(RHess,qp_f,qp_A,qp_B,[],[],qp_XN,0,-1,Nlconst, ...
               size(qp_A,1),size(Z,2),OPTIONS,sqpopt,ACTIND_lu);

    % Second order correction step
    SD_c = Y*py_c+Z*pz_c;
else
    % Compute the correction step in null space, pz_c
    qp_f = Z'*gf+2*pdamp*w+RHess*pz;
    qp_A = [-Z;Z];
    qp_B = [-lb+XOUT+2*Ypy+Zpz;ub-XOUT-2*Ypy-Zpz];
    qp_XN = sparse(size(Z,2),1);
    [pz_c,lambda_lu_c,exitflagqp,outputqp,howqp,ACTIND_lu_c] = ...
        qpsub2(RHess,qp_f,qp_A,qp_B,[],[],qp_XN,0,-1,Nlconst, ...
               size(qp_A,1),size(Z,2),OPTIONS,sqpopt,ACTIND_lu);

    % Second order correction step
    SD_c = Ypy+Z*pz_c;
end

if changeflag
    x = IPMULT*(XOUT+Ypy+Zpz+SD_c);
else
    x = XOUT+Ypy+Zpz+SD_c;
end
f = feval(funfcn{3},x,varargin{:});
if constflag
    [nctmp,nceqtmp] = feval(confcn{3},x,varargin{:});
else
    nceqtmp = [];
end
c = [Aeq*x-Beq; nceqtmp(:); A*x-B];
ga = [abs(c(1:eq)); c(eq+1:end)];
if ~isempty(c)
    mg = max(max(ga),0);
else
    mg = 0;
end

temp = ga;
temp(isinf(ga) | (ga < 0)) = 0;
MERIT = f+sum(OLDLAMBDA.*temp);
if isempty(infeas)
    DGMATL = param+gf'*SD_c;
    MeritCondition = (MERIT > MATL+0.1*DGMATL);
else
    if mg > 0
        MERIT2 = mg;
    elseif f >=0 
        MERIT2 = -1/(f+1);
    else 
        MERIT2 = 0;
    end
    if ~infeas & f < 0
        MERIT2 = MERIT2 + f - 1;
    end
    MeritCondition = ((MERIT2 > param) & (MERIT > MATL));
end

if MeritCondition
    SOC_flag = -1; % correction failed
else
    SOC_flag = 1; % correction succeeded
    ACTIND_lu = ACTIND_lu_c;
    disp('Second Order Correction Succeeded')
end

