%%-------------------- Experiment III ---------------------------------%%
%% From the article 'Large-Scale Quasi-Newton Trust-Region Methods
%% With Low Dimensional Linear Equality Constraints' J.J. Brust, R.F. Marcia,
%% C.G. Petra

% Comparison of solvers on the Rosenbrock function, subject to linear
% equality constraints from netlib problems: http://www.netlib.org/lp/data/

% Initial version: J.B., 07/03/18

%{
    Change log:

    07/21/17, J.B., Computes solutions with proposed solver
    "LMTR_EIG_inf_2_DENSE_G1_CONST_V2" and Matlab's large-scale
    trust-region optimization solver "fmincon".

    08/24/17, J.B., Comparison of two versions of the method. [LMTR_..._V2]
    uses a shape-changing norm on g, whereas [LMTR_..._V3] applies the same
    norm to an oblique projection of g.

    08/28/17, J.B., Computations using a modified solver based on a
    projected gradient stopping condition [LMTR_..._V4].

    08/28/17, J.B., Computations using a modified solver based on 
    approximations to the Lagrangian [LMTR_..._V5].

    09/07/17, J.B., Test on the l-2 based constrained solver

    10/17/17, J.B, Compact representation and eigenvalues of B^{-1}W solver

    02/05/18, J.B, Minor edits to prepare numerical results for manuscript.

    06/05/18, J.B., Use modifed solver that computes performance ratio
    'rho' to accept steps. Tighter convergence bounds, too. 

    07/02/18, J.B., This version allows for comparison for large-scale n.

    07/03/18, J.B., Test on the rosen objective with linear
    equality constaints selected from netlib.
 
    07/04/18, J.B., Comparison with other solvers. RSQP, fmincon-interior,
    possibly fmincon-trust-region

    07/19/18, J.B., Preparation for release
%}


clc;
clear;

addpath ../main
addpath ../auxiliary
addpath ../netlib/readmps

wtest       = warning('off','all');
currentpath = pwd;

datapath    = fullfile(currentpath,'..','/data/');
figpath     = fullfile(currentpath,'..','/figs/');
probpath    = fullfile(currentpath,'..','/netlib/');

rng(090317);

fprintf('---------------- EXPERIMENT III ------------------------\n');

%%-------------------- Netlib problems ----------------------------------%
problems = {'fit1d.mps',...
            'fit2d.mps',...
            'd6cube.mps',...
            'scsd1.mps',...
            'scsd6.mps',...
            'scsd8.mps'};

numn      = length(problems);
numsol    = 3;
                
%%------------------- Generate data stores ---------------------------- %%

exs             = zeros(numn,numsol);
convs           = zeros(numn,numsol);
nbs             = zeros(numn,numsol);
its             = zeros(numn,numsol);
times           = zeros(numn,numsol);
numf            = zeros(numn,numsol);
numg            = zeros(numn,numsol);

rks             = zeros(numn,1);

%%------------------ Solver paramters --------------------------------- %%
                                                                
%% L2-Const, SC-Const
options_const.storedat  = 0;
options_const.hastrrad  = 1;
options_const.ctol      = 1e-5;
options_const.btol      = 1e-10;
options_const.dflag     = 0;
options_const.gtol      = 1e-5;

%% fmincon
% interior-point algorithm
options_fmin_int        = optimoptions('fmincon',...   
                                        'GradObj',      'on',...
                                        'Display',      'off',...
                                        'Algorithm',    'interior-point', ...
                                        'Hessian',      'lbfgs',...
                                        'MaxIter',      1e6, ...
                                        'MaxFunEvals',  1e6, ...
                                        'TolX',         1e-10);
% trust-region-reflective algorithm
options_fmin_tr        = optimoptions('fmincon',...   
                                        'GradObj',      'on',...
                                        'Display',      'off',...
                                        'Algorithm',    'trust-region-reflective', ...                                        
                                        'MaxIter',      1e6, ...
                                        'MaxFunEvals',  1e6);
                                        
% Convergence tolerance: If nomr(Proj(g)) < ctol and norm(Ax-b) < btol,
% then converged.
ctol                    = 1.e-5;                                    
                                    
%%----------------------- Objective function --------------------------- %%
fun             = @(x)( rosen(x) );   

fprintf('\n**********************\nRunning Solver comparison\n');
%fprintf('n\t time RSQP\t time SC\t time L2\t time  fmincon-I\t fmincon-TR\n');
fprintf('n\t m\t rnk\t t-TRSC\t t-TRL2\t t-fmin\n');
for i = 1:numn
    
    name    = problems{i};
    prob    = readmps(fullfile(probpath,name));
    
    A       = prob.A;
    [mn,n]  = size(A);
    
    rnk     = rank(full(A));                      
    
    b0      = randn(n,1);
    b       = A*b0;
    x0      = A'*((A*A')\b);

    const   = @(x)( const_quad_arg(x,A,b));
    
    btol    = norm(A*x0-b)*10;
    
    fprintf('%i\t %i\t %i\t %s\t %s\t %s\n',n,mn,rnk,'------','------','------');
    
    %--------------------- Parameters (based on problem data) ---------%%
    
    %% L2-Const, SC-Const    
    options_const.btol          = btol;
    trradb                      = norm(x0);
    options_const.trradb        = trradb;   
    options_const_l2            = options_const;
    options_const_l2.maxitroot  = 10;
    options_const_l2.epsroot    = 1e-5;
    
    %% fmincon
    options_fmin_int    = optimoptions(options_fmin_int,...
                                        'TolCon', btol);
    options_fmin_tr     = optimoptions(options_fmin_tr,...
                                        'TolCon', btol);
    
    sidx                = 0;
    %%-------------------- Solver calls --------------------------------%%
    %% RSQP
%     sidx                = sidx + 1;
%     
%     tic;
%     [x,f_rsqp,ex_rsqp,out_rsqp] = ...
%     rsqp(fun,x0,A,b,A,b,[],[],[],options_rsqp);
%     time_rsqp = toc;
%     
%     [f,g]           = fun(x);
%     ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
%     err             = ngp/norm(x);    
%     nb              = norm(A*x-b);
%     
%     exs(i,sidx)     = (err < 1e-5) && (nb < btol);
%     convs(i,sidx)   = err;
%     nbs(i,sidx)     = nb;
%     its(i,sidx)     = out_rsqp.iterations;
%     times(i,sidx)   = time_rsqp;
%     numf(i,sidx)    = out_rsqp.funcCount;
%     numg(i,sidx)    = NaN;
    
    %% SC-CONST
    sidx            = sidx + 1;
    [x,~,outinfo]   = LTRSC_LEC_V1(fun,const,x0,options_const); % V6
    
    [~,g]           = fun(x);
    ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
    err             = ngp/norm(x);    
    nb              = norm(A*x-b);
    
    exs(i,sidx)     = (err < ctol) && (nb < btol);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = outinfo.numit;
    times(i,sidx)   = outinfo.tcpu;
    numf(i,sidx)    = outinfo.numg;
    numg(i,sidx)    = outinfo.numf;
    
    %% L2-CONST
    sidx            = sidx + 1;
    [x,~,outinfo]   = LTRL2_LEC_V1(fun,const,x0,options_const_l2);    % _V1
    
    [~,g]           = fun(x);
    ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
    err             = ngp/norm(x);    
    nb              = norm(A*x-b);
    
    exs(i,sidx)     = (err < ctol) && (nb < btol);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = outinfo.numit;
    times(i,sidx)   = outinfo.tcpu;
    numf(i,sidx)    = outinfo.numg;
    numg(i,sidx)    = outinfo.numf;
    
    %% fmincon-I
    sidx            = sidx + 1;
    tic;
    [x,~,ex_fmin,out_fmin] = fmincon(fun,x0,[],[],A,b,...
        [],[],[],options_fmin_int);
    time_fmincon = toc;
    
    [f,g]           = fun(x);
    ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
    err             = ngp/norm(x);    
    nb              = norm(A*x-b);
    
    exs(i,sidx)     = (err < ctol) && (nb < btol);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = out_fmin.iterations;
    times(i,sidx)   = time_fmincon;
    numf(i,sidx)    = out_fmin.funcCount;
    numg(i,sidx)    = out_fmin.funcCount;
    
    fprintf('%i\t %i\t %i\t %6.4f\t %6.4f\t %6.4f\n',n,mn,rnk,times(i,1),times(i,2),times(i,3));
    
%     %% fmincon-TR
%     sidx            = sidx + 1;
%     tic;
%     [x,f,ex_fmin,out_fmin] = fmincon(fun,x0,[],[],A,b,...
%         [],[],[],options_fmin_tr);
%     time_fmincon = toc;
%     
%     [f,g]           = fun(x);
%     ngp             = sqrt(abs(g'*g -g'*A'*((A*A')\(A*g))));
%     err             = ngp/norm(x);    
%     nb              = norm(A*x-b);
%     
%     exs(i,sidx)     = (err < 1e-5) && (nb < btol);
%     convs(i,sidx)   = err;
%     nbs(i,sidx)     = nb;
%     its(i,sidx)     = out_fmin.iterations;
%     times(i,sidx)   = time_fmincon;
%     numf(i,sidx)    = out_fmin.funcCount;
%     numg(i,sidx)    = out_fmin.funcCount;
%     
%     %%-------------- Display times -------------------------------------%%
%     fprintf('%d\t %.3e\t %.3e\t %.3e\t %.3e \t %.3e\n',n,...
%         times(i,1),...
%         times(i,2),...
%         times(i,3),...
%         times(i,4),...
%         times(i,5));    
end

%%------------------- Comparison plots ---------------------------------%%

%   'RSQP',...
leg={'TR-$(\mathbf{P},\infty)$',...
        'TR-$\ell_2$',...       
        'fmincon-I-ldl'};%,...
        %'FMINCON-TR'
                
types.colors    = ['b' 'r' 'y' ]; %'k' 'y'
types.lines     = {'-', '-.', '-', }; %'-',   '-'
types.markers   = ['o' 'o' 'o' ]; %'s' 's'
                
indAlg          = [1 2 3 ]; %4 5

perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_III.eps'));


perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_III.eps'));

save(fullfile(datapath,'EXPERIMENT_III'),'exs','convs','nbs','its','times','numf','numg');

close ALL;

% Restore warning settings
warning(wtest);

