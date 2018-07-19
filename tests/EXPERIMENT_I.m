%% EXPERIMENT I from the article 'Large-Scale Quasi-Newton Trust-Region Methods
%% With Low Dimensional Linear Equality Constraints' J.J. Brust, R.F. Marcia,
%% C.G. Petra
%{
    Here the test is on quadratic programming problems

    min 1/2 x' Q x + c' x,  subject to,     Ax = b,

    where Q p.s.d (nxn) and n varies.

    Initial version: 07/16/18, J.B.

    07/18/17, J.B., Save experiment data, and plots to folders 'data' and
    'figs'.
%}

clc;
clear;

addpath ../main
addpath ../auxiliary
addpath ../solvers/rsqp
addpath ../solvers/pdco

wtest = warning('off','all');

rng(090317);

fprintf('---------------- EXPERIMENT I ------------------------\n');

% Branch by 'small' or 'larger'
islarge         = input(['Please choose test size: Small=0  Large=1.\n', ...
                            'Then hit <Enter> :']);
                           
if ~islarge
    ns = 50:50:500;
else
    ns = 1000:100:2000;
end

%------------------------- Data storage ---------------------------------
mm              = 10;
numn            = length(ns);
numsol          = 6;

exs             = zeros(numn,numsol);
convs           = zeros(numn,numsol);
nbs             = zeros(numn,numsol);
its             = zeros(numn,numsol);
times           = zeros(numn,numsol);
numf            = zeros(numn,numsol);
numg            = zeros(numn,numsol);

idxsol          = 0;

ctol1           = 1e-6; % Convergence: abs((f^* - f^k)/f^*) < ctol1
ctol2           = 1e-9; % Feasibility: norm(Ax^k-b) < ctol2

%------------------------- Solver parameters ------------------------------
%% Solver 1: PDCO
options_pdco            = pdcoSet;
options_pdco.mu0        = 0e-0;  % An absolute value
options_pdco.Method     = 21;    % (1=chol  2=QR  3=LSMR  4=MINRES  21=SQD(LU)  22=SQD(MA57))
options_pdco.LSQRatol1  = 1e-8;  % For LPs, LSQR must solve quite accurately
options_pdco.wait       = 0;     % Allow options to be reviewed before solve
options_pdco.Print      = 0;

%% Solver 2: RSQP 
options_rsqp            = optimset2('GradObj','on', ...
                                    'Display','off',...
                                    'MaxIter', 1e6);

%% Solver(s) 3,4: LMTR_SC, LMTR_L2
options_sc.storedat = 0;
options_sc.btol     = 1e-10;
options_sc.dflag    = 0;
options_sc.gtol     = 1e-5;

%% Solver 5: Fmincon Interior Ldl
options_fmin_ldl    = optimoptions('fmincon',...   
                                            'GradObj',      'on',...
                                            'Display',      'off',...
                                            'Algorithm',    'interior-point', ...
                                            'Hessian',      'lbfgs',...
                                            'MaxIter',      1e6, ...
                                            'MaxFunEvals',  1e6, ...
                                            'TolX',         1e-10,...
                                            'SubproblemAlgorithm','ldl-factorization'); %% cg
                                        
%% Solver 6: Fmincon Interior Cg
options_fmin_cg    = optimoptions('fmincon',...   
                                            'GradObj',      'on',...
                                            'Display',      'off',...
                                            'Algorithm',    'interior-point', ...
                                            'Hessian',      'lbfgs',...
                                            'MaxIter',      1e6, ...
                                            'MaxFunEvals',  1e6, ...
                                            'TolX',         1e-10,...
                                            'SubproblemAlgorithm','cg');                                        

fprintf('\n**********************\nRunning Solver comparison\n');
fprintf('n\t time PDCO\t time RSQP\t time LTRSC\t time LTRL2 \t time fmin-ldl \t time fmin-cg\n');

for i = 1:numn

    % Problem data
    n   = ns(i);

    Q1  = randn(n,n); 
    Q   = Q1'*Q1;
    A   = randn(mm,n);  
    b0  = randn(n,1);    
    b   = A*b0;
    c   = randn(n,1);
    x0  = A'*((A*A')\b);
    
    fun         = @(x)( object_quad_arg(x,c,Q));
    const_sc    = @(x)( const_quad_arg(x,A,b));
    
    % --------------------- Parameters (based on problem data) -------------
    %% Solver 1: PDCO
    em      = ones(mm,1);
    en      = ones(n,1);
    zn      = zeros(n,1);
    bigbnd  = 1e+30;
    gamma   = 1e-5;              % Primal regularization
    delta   = 1e-5;              % 1e-3 or 1e-4 for LP;  1 for Least squares.
    d1      = gamma;             % Can be scalar if D1 = d1*I.
    d2      = delta*em;
    bl_pdco = -Inf.*en;
    bu_pdco = Inf.*en;
    y0      = zeros(mm,1);        % Initial y
    z0      = ones(n,1)/n;       % Initial z
    xsize   = 5;                 % Estimate of norm(x,inf) at solution
    zsize   = 8;                 % Estimate of norm(z,inf) at solution
    
    sQ       = sparse(Q);
    obj_pdco = @(x)(deal(0.5*(x'*Q*x) + c'*x, Q*x + c, sQ));
    
    % ----------------------- Analytic solution ---------------------------
    Qi      = Q\[A' c];
    Qc      = Qi(:,end);
    QAt     = Qi(:,1:(end-1));
    xopt    = -Qc + QAt*((A*(QAt))\(b+A*(Qc)));
    fopt    = fun(xopt);
    
    %% Solver(s) 4,5: LMTR_SC, LMTR_L2
    options_sc.trradb      = norm(x0);
    
    options_l2            = options_sc;
    options_l2.maxitroot  = 10;
    options_l2.epsroot    = 1e-5;
    
     sidx                = 0;
     
    % ---------------------- Solver calls ---------------------------------
    %% Solver 1
    sidx                = sidx + 1;
    [x_pdco,y_pdco,z_pdco,inform_pdco,PDitns_pdco,CGitns_pdco,time_pdco] = ...
    pdco(obj_pdco,A,b,bl_pdco,bu_pdco,d1,d2,options_pdco,x0,y0,z0,xsize,zsize );
    
    err             = abs(fopt-fun(x_pdco))/abs(fopt);
    nb              = norm(A*x_pdco-b);
    
    exs(i,sidx)     = (err < ctol1) && (nb < ctol2);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = PDitns_pdco;
    times(i,sidx)   = time_pdco;
    numf(i,sidx)    = NaN;
    numg(i,sidx)    = NaN;

    %% Solver 2
    sidx            = sidx + 1;
    tic;
    [x_rsqp,f_rsqp,ex_rsqp,out_rsqp] = ...
    rsqp(fun,x0,A,b,A,b,[],[],[],options_rsqp);
    time_rsqp = toc;
    
    err             = abs(fopt-f_rsqp)/abs(fopt);
    nb              = norm(A*x_rsqp-b);
    
    exs(i,sidx)     = (err < ctol1) && (nb < ctol2);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = out_rsqp.iterations;
    times(i,sidx)   = time_rsqp;
    numf(i,sidx)    = out_rsqp.funcCount;
    numg(i,sidx)    = NaN;
    
    %% Solver 3
    sidx    = sidx + 1;
    [x_sc,f_sc,out_sc] = LTRSC_LEC_V1(fun,const_sc,x0,options_sc); % LTRSC_LEC
    
    err             = abs(fopt-f_sc)/abs(fopt);
    nb              = norm(A*x_sc-b);
    
    exs(i,sidx)     = (err < ctol1) && (nb < ctol2);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = out_sc.numit;
    times(i,sidx)   = out_sc.tcpu;
    numf(i,sidx)    = out_sc.numf;
    numg(i,sidx)    = out_sc.numg;
     
    %% Solver 4
    sidx    = sidx + 1;
    [x_l2,f_l2,out_l2] = LTRL2_LEC_V1(fun,const_sc,x0,options_l2); % LTRL2_LEC
    
    err             = abs(fopt-f_l2)/abs(fopt);
    nb              = norm(A*x_l2-b);
    
    exs(i,sidx)     = (err < ctol1) && (nb < ctol2);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = out_l2.numit;
    times(i,sidx)   = out_l2.tcpu;
    numf(i,sidx)    = out_l2.numf;
    numg(i,sidx)    = out_l2.numg;
    
    %% Solver 5: Fmincon interior-ldl 
    sidx    = sidx + 1;
    tic;
    [x_fmin,f_fmin,ex_fmin,out_fmin] = fmincon(fun,x0,[],[],A,b,...
        [],[],[],options_fmin_ldl);
    time_fmincon_ldl = toc;
    
    err              = abs(fopt-f_fmin)/abs(fopt);
    nb               = norm(A*x_fmin-b);
    
    exs(i,sidx)     = (err < ctol1) && (nb < ctol2);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = out_fmin.iterations;
    times(i,sidx)   = time_fmincon_ldl;
    numf(i,sidx)    = out_fmin.funcCount;
    numg(i,sidx)    = out_fmin.funcCount;
    
    %% Solver 6: Fmincon interior-cg
    sidx    = sidx + 1;
    tic;
    [x_fmin,f_fmin,ex_fmin,out_fmin] = fmincon(fun,x0,[],[],A,b,...
        [],[],[],options_fmin_cg);
    time_fmincon_cg  = toc;
    
    err              = abs(fopt-f_fmin)/abs(fopt);
    nb               = norm(A*x_fmin-b);
    
    exs(i,sidx)     = (err < ctol1) && (nb < ctol2);
    convs(i,sidx)   = err;
    nbs(i,sidx)     = nb;
    its(i,sidx)     = out_fmin.iterations;
    times(i,sidx)   = time_fmincon_cg;
    numf(i,sidx)    = out_fmin.funcCount;
    numg(i,sidx)    = out_fmin.funcCount;

    fprintf('%d\t %.3e\t %.3e\t %.3e\t %.3e \t %.3e \t %.3e\n',n,time_pdco,time_rsqp,...
        out_sc.tcpu,out_l2.tcpu,time_fmincon_ldl,time_fmincon_cg);
          
end
 % ----------------------- Outputs ---------------------------
leg={   'PDCO',...
        'RSQP',...
        'TR-$(\mathbf{P},\infty)$',...
        'TR-$\ell_2$',...       
        'fmincon-I-ldl', ...
        'fmincon-I-cg'};
                
types.colors    = ['k' 'm' 'b' 'r' 'y' 'g'];
types.lines     = {'-', '-.', '-', '-.', '-', '-.'};
types.markers   = ['o' 'o' 'o' 'o' 'o' 'o'];
                
indAlg          = [1 2 3 4 5 6];

currentpath     = pwd;
datapath        = fullfile(currentpath,'..','/data/');
figpath         = fullfile(currentpath,'..','/figs/');

if islarge == true
    
    perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types);    
    print(gcf, '-dpsc2', fullfile(figpath,'time_EX_I_L.eps'));
    
    perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types);
    print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_I_L.eps'));
    
    perf(exs(:,indAlg),numf(:,indAlg),leg(indAlg),1,types);
    print(gcf, '-dpsc2', fullfile(figpath,'funcs_EX_I_L.eps'));
    
    save(fullfile(datapath,'EXPERIMENT_I_L'),'exs','convs','nbs','its','times','numf','numg');    
else
    
    perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types);    
    print(gcf, '-dpsc2', fullfile(figpath,'time_EX_I_S.eps'));
    
    perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types);
    print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_I_S.eps'));
    
    perf(exs(:,indAlg),numf(:,indAlg),leg(indAlg),1,types);
    print(gcf, '-dpsc2', fullfile(figpath,'funcs_EX_I_S.eps'));
    
    save(fullfile(datapath,'EXPERIMENT_I_S'),'exs','convs','nbs','its','times','numf','numg');
end

close ALL;
 
warning(wtest);

