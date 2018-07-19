%%---------------------------- EXPERIMENT II----------------------------%%
%% From the article 'Large-Scale Quasi-Newton Trust-Region Methods
%% With Low Dimensional Linear Equality Constraints' J.J. Brust, R.F. Marcia,
%% C.G. Petra
%{
    Here the test is on quadratic programming problems

    min 1/2 x' H x + c' x,  subject to,     Ax = b,

    where H p.s.d (nxn), with large values of n. Instead of computing a
    dense Hessian of the objective function it is defined as

    H = phi.I + Q Q^T,
    
    where Q in R^{n x r} and phi > 0.
    
    Initial version: 07/16/18, J.B.

    07/17/18, J.B., Application of trust-region methods that use the
    performance metric rho.

    07/18/18, J.B., Storing outputs in separate folders.

    ----------------------------------------------------------------------
    NOTE: For possibly faster methods use LTRSC_LEC, and LTRL2_LEC
%}

clc;
clear;

addpath ../main
addpath ../auxiliary

wtest = warning('off','all');

currentpath     = pwd;
datapath        = fullfile(currentpath,'..','/data/');
figpath         = fullfile(currentpath,'..','/figs/');

rng(090317);

fprintf('---------------- EXPERIMENT II ------------------------\n');

%------------------------- Data storage ---------------------------------
ns              = [1e4;1e5;1e6;1e7];

numn            = length(ns);
mm              = 10;

% Only 4 solvers are compared, because PDCO, RSQP take exceedingly long
% on large-scale problems.
numsol          = 4;

%%------------------- Generate data stores ---------------------------- %%
exs             = zeros(numn,numsol);
convs           = zeros(numn,numsol);
nbs             = zeros(numn,numsol);
its             = zeros(numn,numsol);
times           = zeros(numn,numsol);
numf            = zeros(numn,numsol);
numg            = zeros(numn,numsol);

%%------------------ Solver paramters --------------------------------- %%

%% L2-Const, SC-Const
options_const.storedat  = 0;
options_const.gtol      = 5e-5;
options_const.btol      = 5e-9;
options_const.dflag     = 0;

%% fmincon
% interior-point algorithm
options_fmin_ldl        = optimoptions('fmincon',...
    'GradObj',      'on',...
    'Display',      'off',...
    'Algorithm',    'interior-point', ...
    'Hessian',      'lbfgs',...
    'MaxIter',      500, ... % 1e6
    'MaxFunEvals',  1e6, ...
    'TolX',         1e-10);
% Interior-Cg out of memory
options_fmin_cg        = optimoptions('fmincon',...
    'GradObj',      'on',...
    'Display',      'off',...
    'Algorithm',    'interior-point', ...
    'SubproblemAlgorithm','cg',...
    'MaxIter',      500, ... % 1e6
    'MaxFunEvals',  1e6);


fprintf('\n**********************\nRunning Large-Scale Comparison\n');
fprintf('n\t time LTRSC\t time LTRL2 \t time fmin-ldl \t fmin-cg\n');

for i = 1:numn
    
    n   = ns(i);
    
    phi     = abs(randn(1));
    phii    = 1/phi;
    r       = 5;
    Ir      = eye(r);
    
    Q   = randn(n,r);
    
    A   = randn(mm,n);
    x0  = randn(n,1);
    b   = A*x0;
    c   = randn(n,1);
    
    fun     = @(x)( object_quad_arg_L(phi,x,c,Q));
    
    Qi      = phii.*[A' c] - Q*((phi*phi.*Ir + phi.*Q'*Q)\(Q'*[A' c]));
    Qc      = Qi(:,end);
    QAt     = Qi(:,1:(end-1));
    xopt    = -Qc + QAt*((A*(QAt))\(b+A*(Qc)));
    fopt    = fun(xopt);
    
    const   = @(x)( const_quad_arg(x,A,b));
    
    trradb  = norm(x0);
    
    options_const.trradb   = trradb;
    
    options_const_l2            = options_const;
    options_const_l2.maxitroot  = 10;
    options_const_l2.epsroot    = 1e-5;
        
    sidx                = 0;
    
    %% TR-SC
    sidx                  = sidx + 1;
    [x,f,outinfo]         = LTRSC_LEC_V1(fun,const,x0,options_const); % LTRSC_LEC
    
    err                  = abs(fopt-f)/abs(fopt);
    nb                   = norm(A*x-b);
    
    exs(i,sidx)          = (err < 1e-7) && (nb < 1e-8);
    convs(i,sidx)        = err;
    nbs(i,sidx)          = nb;
    its(i,sidx)          = outinfo.numit;
    times(i,sidx)        = outinfo.tcpu;
    numf(i,sidx)         = outinfo.numf;
    numg(i,sidx)         = outinfo.numg;
    fprintf('%d\t %.3e\t ------ \t ------ \t ------\n',n,times(i,1));
    
    %% TR-L2
    sidx                    = sidx + 1;
    [x,f,outinfo]           = LTRL2_LEC_V1(fun,const,x0,options_const_l2);  % LTRL2_LEC
    
    err                     = abs(fopt-f)/abs(fopt);
    nb                      = norm(A*x-b);
    
    exs(i,sidx)             = (err < 1e-7) && (nb < 1e-8);
    convs(i,sidx)           = err;
    nbs(i,sidx)             = nb;
    its(i,sidx)             = outinfo.numit;
    times(i,sidx)           = outinfo.tcpu;
    numf(i,sidx)            = outinfo.numf;
    numg(i,sidx)            = outinfo.numg;
    
    fprintf('%d\t %.3e\t %.3e \t ------ \t ------\n',n,times(i,1),times(i,2));
    
    exs(i,sidx+1)           = -1;
    exs(i,sidx+2)           = -1;
    
    if n <= 1e5
        %% Fmincon-LDL
        sidx            = sidx + 1;
        try
            tic;
            [x,f,ex_fmin,out_fmin] = fmincon(fun,x0,[],[],A,b,...
                [],[],[],options_fmin_ldl);
            time_fmincon = toc;
            
            err                     = abs(fopt-f)/abs(fopt);
            nb                      = norm(A*x-b);
            
            exs(i,sidx)             = (err < 1e-7) && (nb < 1e-8);
            convs(i,sidx)           = err;
            nbs(i,sidx)             = nb;
            its(i,sidx)             = out_fmin.iterations;
            times(i,sidx)           = time_fmincon;
            numf(i,sidx)            = out_fmin.funcCount;
            numg(i,sidx)            = out_fmin.funcCount;
            
            fprintf('%d\t %.3e\t %.3e \t %.3e \t ------ \n',n,times(i,1),times(i,2),times(i,3));
        catch
            exs(i,sidx)         = -1;
        end
        
        %% Fmincon-CG
        sidx            = sidx + 1;
        try
            tic;
            [x,f,ex_fmin,out_fmin] = fmincon(fun,x0,[],[],A,b,...
                [],[],[],options_fmin_cg);
            time_fmincon = toc;
            
            err                     = abs(fopt-f)/abs(fopt);
            nb                      = norm(A*x-b);
            
            exs(i,sidx)             = (err < 1e-7) && (nb < 1e-8);
            convs(i,sidx)           = err;
            nbs(i,sidx)             = nb;
            its(i,sidx)             = out_fmin.iterations;
            times(i,sidx)           = time_fmincon;
            numf(i,sidx)            = out_fmin.funcCount;
            numg(i,sidx)            = out_fmin.funcCount;
            
            fprintf('%d\t %.3e\t %.3e \t %.3e \t %.3e \n',n,times(i,1),times(i,2),times(i,3),...
                times(i,4));
            
        catch
            exs(i,sidx)         = -1;
        end        
    end    
end

leg={   'TR-$(\mathbf{P},\infty)$',...
        'TR-$\ell_2$',...
        'fmincon-I-ldl', ...
        'fmincon-I-cg'};

types.colors    = ['b' 'r' 'y' 'g'];
types.lines     = {'-', '-.', '-', '-.'};
types.markers   = ['o' 'o' 'o' 'o'];

indAlg          = [1 2 3 4];
types.colors    = types.colors(indAlg);

% Modify location of legend to have the same output as in article, i.e.
% change loaction to 'NorthEast' from 'SouthEast'.
perf(exs(:,indAlg),times(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'time_EX_II.eps'));


perf(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types);
print(gcf, '-dpsc2', fullfile(figpath,'iter_EX_II.eps'));

save(fullfile(datapath,'EXPERIMENT_II'),'exs','convs','nbs','its','times','numf','numg');

close ALL;

% Restore warning settings
warning(wtest);
