function [x,f,outinfo] = LTRL2_LEC(fun,const,x0,params)
% LTRL2_LEC - Limited-Memory Trust-Region Shape-Changing Norm for
% Linear Equality Constrained Problems
%
% For details, see: 
% Large-Scale Quasi-Newton Trust-Region Methods With Low Dimensional
% Linear Equality Constraints, J.J. Brust, R.F. Marcia, C.G. Petra, 2018
%
% Method developed for problems  
%
%       min f(x), s.t Ax = b_h,
%
% where
%
% f := fun(x), Ax - b_h := const(x).
%
% Solution approach is a sequential quadratic programming (SQP) technique,
% where at a given iterate
%
%           x_{k+1} = x_k + s_k, 
%
% the step s_k is computed by
%
% s_k = argmin {Q(s) := 1/2 s' B s + s'g}, subject to:
%
%    As = 0, || s ||_{P,infty} <= Delta_k,
%
% where g = f'(x_k), B approx f''(x_k), A = h'(x_k), b = A(x_k)-b_h approx 0. 
% The solution to the trust-region constraint is computed by solving
% the equation:
%
% phi(sigma) = 1/|| s(sigma) ||_2 - 1/Delta_k.
%
% At a solution it holds phi(sigma*)= 0 and sigma* > 0.
%
% Here Delta_k is the "trust-region" radius and B = B_k is the L-BFGS
% compact representation: 
%
%       B = delta.I + V M V',   and     B^{-1} = H = gamma.I - V N V',
%
% The notation corresponding to the article is V (code) = \hat{\Psi} (article) and 
% delta (code) = gamma (article), \hat{M} (code) = -\hat{M} (article).
% The initial parameters are defined in a struct params. 
% The output information is contained in a struct outinfo.
%
% params contains the following fields {default values}:
%   m       - maximum number of stored vector pairs (s,y) {5}
%   [gtol    - tolerance on L2-norm of gradient ||g||<gtol*max(1,||x||)
%   {1e-5}]
%   ftol    - tolerance for abs(f_k-f_{k+1})
%   maxit   - maximum number of iterartions {100000}
%   ranktol - tolerance threshold for establishing rank of V {1e-7}
%   dflag   - display parameter {1}:
%               1 - display every iteration;
%               0 - no display.
%   storedat- flag whether or not to store data
%               1 - store:
%                       - Delta_k,                       
%                       - f_k
%                       - Convergence criterion,
%   trradb  - trrad_k = trradb + trrad (for linearly constrained
%   problems)
%   btol    - ||h_k - b_h||<btol
%   maxitroot - maximum number of iterations for the root finding problem
%   epsroot - error in the root finding process
%
% outinfo contains the following fields:    
%   ex      - exitflag:
%               1 - Projected gradient satisfies termination criterion
%              -1 - TR radius is too small
%              -2 - exceeds maximum number of iterations
%              >1 - line search failed, see cvsrch.m
%   numit   - number of succesful TR iterations
%   numf    - number of function evaluations
%   numc    - number of constraint evaluations
%   numg    - number of gradient evaluations
%   numcg   - number of constraint gradient evaluations
%   numrst  - number of restarts when V'*s is of low accuracy
%   tcpu    - CPU time of algorithm execution
%   tract   - number of iterations when TR is active
%   trrej   - number of iterations when initial trial step is rejected
%   params  - input paramaters
% 
% See also LMTR_out, svsrch
%
%
% This code is distributed under the terms of the GNU General Public
% License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that (A) this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software, and (B) that the corresponding 
% article is appropriately cited. 
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.

%{

J.B., 20/01/17 
Modification to allow for different choices of gamma_perp in the
trust-region subproblem. This is only applicable if the quasi-Newton step
is rejected. The L-BFGS matrix has form

    B_k = B_0 + V M V',

             | gamma.I          0     |
    B_0 =  P |   0        gamma_perp.I| P'

J.B., 08/23/17
Apply shape-changing norm to 

	||Wg||_{P,inf} 	instead of ||g||_{P,inf},

	where

	W = I - A'(A B^{-1} A')^{-1} A B^{-1} (oblique projection)

J.B., 08/24/17
Modification of stopping criterion. Instead of the true solution, seek a
point for which || Proj(g) || = 0.

J.B., 07/02/18
Modification for rank degenerate problems. This is based on an LDL'
factorization of AA'.

J.B., 07/16/18
Removal of warnings.
%}

% Read input parameters
if nargin<3
    params = struct;
end;
  
if isfield(params,'m') && ~isempty(params.m)
    m = params.m;
else
    m = 5;
end;

if isfield(params,'ftol') && ~isempty(params.ftol)
    ftol = params.ftol;
else
    ftol = 1e-11;
end;

if isfield(params,'gtol') && ~isempty(params.gtol)
    gtol = params.gtol;
else
    gtol = 1e-5;
end;

if isfield(params,'maxit') && ~isempty(params.maxit)
    maxit = params.maxit;
else
    maxit = 100000;
end;

if isfield(params,'ranktol') && ~isempty(params.ranktol)
    ranktol = params.ranktol;
else
    ranktol = 1e-10;
end;

if isfield(params,'dflag') && ~isempty(params.dflag)
    dflag = params.dflag;
else
    dflag = 1;
end;

if isfield(params,'storedat') && ~isempty(params.storedat)
    storedat    = params.storedat;
    trrads      = zeros(maxit,1); % store Delta_k,
    funs        = zeros(maxit,1); % f_k,
    converr     = zeros(maxit,1); % || Proj(g) || / max(x,1) 
else
    storedat    = 0;    
end;

% if isfield(params,'hastrrad') && ~isempty(params.hastrrad)
%     hastrrad = params.hastrrad;
% else
%     hastrrad = 0;
% end;

if isfield(params,'trradb') && ~isempty(params.trradb)
    trradb = params.trradb;
else
    trradb = 0;
end;

% if isfield(params,'fopt') && ~isempty(params.fopt)
%     fopt = params.fopt;
% else
%     fopt = 0;
% end;
% 
% if isfield(params,'ctol') && ~isempty(params.ctol)
%     ctol = params.ctol;
% else
%     ctol = 1e-5;
% end;

if isfield(params,'btol') && ~isempty(params.btol)
    btol = params.btol;
else
    btol = 1e-5;
end;

if isfield(params,'maxitroot') && ~isempty(params.maxitroot)
    maxitroot = params.maxitroot;
else
    maxitroot = 10;
end;

if isfield(params,'epsroot') && ~isempty(params.epsroot)
    epsroot = params.epsroot;
else
    epsroot = 1e-5;
end;
% Set trust region parameters using the following pseudo-code
%
% if ratio>=tau0 
%   trial step is accepted 
% end
% if ratio>=tau1 and ||s*||>=c3*trrad
%   trrad=c4*trrad 
% elseif ratio<tau2
%   trrad=max(c1*||s||,c2*trrad)
% end
% if trrad<trrad_min
%   return
% end
%
% Accept the trial step also if (ft-f)/abs(f)<ftol && (ratio>0)

tau0    = 0;                     
tau1    = 0.25;                  
tau2    = 0.75;                  
c1      = 0.5;                     
c2      = 0.25;                    
c3      = 0.8;
c4      = 2;
trrad_min = 1e-15; 

% Set parameters for More-Thuente line search procedure cvsrch

gtolls  = 0.9;
ftolls  = 1e-4;
xtolls  = 1e-16;
maxfev  = 20;
stpmin  = 1e-20;
stpmax  = 1e20;

% Set linsolve options for solving triangular linear systems
opts1.UT        = true;
opts1.RECT      = false;
opts1.TRANSA    = false;

opts2           = opts1;
opts2.TRANSA    = true; % Transposed system

opts3.SYM       = true; % Symmetric system
% Start measuring CPU time
t0=tic;     

%% Memory Allocation 

n=size(x0,1);
m2=m*2;


% Allocate memory for elements of L-BFGS matrix B = delta*I + V*W*V'
V       = zeros(n,m2);  % V=[S Y]
nV      = zeros(m2,1);  % norms of column vectors in V
Vg      = zeros(m2,1);  % product V'*g for the new iterate
Vg_old  = zeros(m2,1);  % product V'*g for the previous iterate
VV      = zeros(m2,m2);  % Gram matrix V'*V
T       = zeros(m,m);  % upper-triangular part of S'*Y
L       = zeros(m,m);  % lower-triangular part of S'*Y with 0 main diagonal
E       = zeros(m,1);  % E=diag(S'*Y)

alpha   = 0; %#ok<NASGU>

% Initialize indexes and counters
numsy   = 0;  % number of stored couples (s,y)
numsy2  = 0;  % numsy2=numsy*2
maskV   = [];  % V(:,maskV)=[S Y]
rankV   = 0;  %#ok<NASGU> % column rank of V
Nflag   = 0;  % indicator, 1 if quasi-Newton step is accepted, 0 if rejected
tract   = 0;  % number of iterations when TR is active
trrej   = 0;  % number of iterations when initial trial step is rejected
it      = 0;  % number of TR iterations
numf    = 0;  % number of function evaluations
numg    = 0;  % number of gradient evaluations
numrst  = 0;  % number of restarts

numskip = 0; % number of update skips, when s'*y<=0
numgp0  = 0; % number of TR subproblems in which g_perp~0

numc    = 0; % number of constraint evaluations
%numcg   = 0; % number of constraint gradient evaluations 

% Constraints

[A,b]       = const(x0);
numc        = numc + 1;
%numcg       = numcg + 1;
mm          = length(b);
maskA       = 1:mm;
rankA       = length(maskA);

AA          = zeros(rankA,rankA); % AA'

VA          = zeros(m2,mm); % V'A'
RA          = zeros(mm,mm); % RA'*RA = AA', Cholesky or LDL' factor of AA'

Uc          = zeros(m2,mm+1); % [gamma.* N] [ V'A' V'g ]

b1          = zeros(mm,1); % b+...
b2          = zeros(mm,1); % (AA'+...)^{-1} b2

Ag          = zeros(mm,1); % A'g
Ag_old      = zeros(mm,1); % A*g previous iteration
Vgradp      = zeros(m2,1); % V'*W*g

Vs          = zeros(m2,1); % V's
As          = zeros(mm,1); % A s
N12         = zeros(m2,m2);

%{ 

    Here 
    
        [gamma.N]^{-1} = [ (1/alpha -1).S'S,               (1/alpha -1).S'Y- alpha.T;
                               (1/alpha -1).Y'S- alpha.T'      -(gamma.D +
                               Y'Y) ],                                                              
    where

    alpha = gam/gam+sigma and gamma := gamma + sigma.
    
%}

%betac   = 0; % B^{-1} B_1 = P diag(lam d1, beta.I) P'


%% Initial check for optimality

% Evaluate function and gradient in the starting point

[f0, g0]    = fun(x0);
numf        = numf+1;
numg        = numg+1;
ng          = norm(g0);

%% Initialization: line search along normalized steepest descent direction,
%% projected onto the linearized feasible set.
% Assume initial feasible point
    
it  = it+1;
x   = x0;
f   = f0;
g   = g0;

% Row scaled version
sA      = sqrt(sum(A.^2,2));
idxnzr  = find(sA > 1e-10);
dsA     = diag(sA);
% Upper triangular matrix, and rank determination of A.
[Ldla, Ddla, pdla]  = ldl(dsA(idxnzr,(idxnzr))\...
                        ((A((idxnzr),:)*A((idxnzr),:)')/...
                        dsA((idxnzr),(idxnzr))),'vector');
                    
dDdla               = diag(Ddla);

% Compute linearly independent columns.
maskldla            = find(dDdla>ranktol); % ^2
rankA               = length(maskldla);  % column rank of Psi
maskA               = idxnzr(pdla(maskldla));  % index of safely linearly independent columns of Psi

RA(1:rankA,1:rankA)     = diag(sqrt(dDdla(maskldla)))*(Ldla(maskldla,maskldla)'*dsA(maskA,maskA));
AA(1:rankA,1:rankA)     = A(maskA,:)*A(maskA,:)';
Ag(1:rankA)             = A(maskA,:)*g;

b2(1:rankA)             = linsolve(RA(1:rankA,1:rankA),...
                                    linsolve(RA(1:rankA,1:rankA),Ag(1:rankA),opts2),...
                                    opts1);

gc                      = A(maskA,:)'*b2(1:rankA);
nb                      = norm(b);

d_                      = gc-g;
ngp                     = sqrt(ng*ng-Ag(1:rankA)'*b2(1:rankA));
d                       = d_./ngp;

conv                    = ngp/max(1,norm(x0));

if dflag==1
    fprintf('\n**********************\nRunning LTRL2-LEC\n');
    fprintf('it\t obj\t\t norm(b)\t norm(Proj(g))) \t norm(dx)\t trrad\n');
    fprintf('%d\t %.3e\t %.3e\t %.3e\t---\t\t ---\n',0,f0,nb,ngp);
end

if ((conv<gtol) && (nb < btol)) % ng<max(1,norm(x0))*gtol
    ex=1;
    x=x0;
    f=f0;
    tcpu=toc(t0);   
    
    outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,numskip,numgp0);
    return;
end

% Shift x, so that it becomes feasible.
if nb > btol
    
    xl      = x -A(maskA,:)'*linsolve(RA(1:rankA,1:rankA),...
                                    linsolve(RA(1:rankA,1:rankA),Ag(1:rankA),opts2),...
                                    opts1);   
else
    
    xl      = x;
    
end  

ns      = 1;
xt      = xl+ns*d;    
ft      = fun(xt);
numf    = numf+1;

% Backtracking
if ft<f  % doubling step length while improvement    
    f=ft;
    ns=ns*2;    
    xt=xl+ns*d;    
    ft=fun(xt);
    numf=numf+1;
    while ft<f
        f=ft;
        ns=ns*2;
        xt=xl+ns*d;
        ft=fun(xt);
        numf=numf+1;
    end
    ns=ns/2;    
else  % halving step length until improvement        
    while ft>=f
        ns=ns/2;        
        if ns<trrad_min
            tcpu=toc(t0);
            ex=-1;
            
            outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,numskip,numgp0);
            return;
        end
        xt=xl+ns*d;
        ft=fun(xt);
        numf=numf+1;
    end
    f=ft;       
end  % line search

g_old           = g;
Ag_old(1:rankA) = Ag(1:rankA);
s               = ns*d;
ns              = norm(s);
%x=x+s;
x               = xl+s;

[~,g]           = fun(x);
numg            = numg+1;
ng              = norm(g);

Ag(1:rankA)     = A(maskA,:)*g;
b2(1:rankA)     = linsolve(RA(1:rankA,1:rankA),...
                                    linsolve(RA(1:rankA,1:rankA),Ag(1:rankA),opts2),...
                                    opts1);

ngp             = sqrt(ng*ng-Ag(1:rankA)'*b2(1:rankA));
conv            = ngp/max(1,norm(x));

[~,b]           = const(x);
numc            = numc + 1;

if (conv>gtol) || (nb > btol)  % norm of gradient is too large, ng>gtol*max(1,norm(x))
    mflag   = 1;  % new iteration is to be made
    y       = g-g_old;
    ny      = norm(y);
    sy      = s'*y;
    
    %{
    if (sy<=1.0e-8*ns*ny)
        error("Initial curvature condition, s'y > 0, violated.");
    end
    %}
    
    % Try More-Thuente line search if positive curvature condition fails
    if (sy<=1.0e-8*ns*ny)
        [x,f,g,ns,exls,numfls] = ...
            cvsrch(fun,n,xl,f0,g_old,d,1,ftolls,gtolls,xtolls,stpmin,stpmax,maxfev); %xl for x0
        
        numf = numf + numfls;
        numg = numg + numfls;        
        
        if (exls>1)  % line search failed
            
            ex      = exls;
            tcpu    = toc(t0);                   
            
            outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,numskip,numgp0);
            return;
        end;            

        s   = x-x0;
        y   = g-g_old;
        ny  = norm(y);
        sy  = s'*y;
        
    end
    
    
else
    mflag = 0;  % problem is solved
end

trrad=2*ns;  % take TR radius as the doubled initial step length
if ns~=1
    tract=tract+1;
end

% Display information about the last iteration
if (dflag==1)
    fprintf('%d\t %.3e\t %.3e\t %.3e\t %.3e\t ---\n',it,f,nb,ngp,ns,trrad+trradb);
end  

if (storedat ~= 0) 
    trrads(it)      = trrad+trradb; % store Delta_k,
    funs(it)        = f; % f_k,
    converr(it)     = conv;
end

%% Main loop
while (mflag==1)    
    %% L-BFGS update    
    if (sy>1.0e-8*ns*ny)  % positive curvature condition holds, update
        delta=ny^2/sy;        
        if numsy<m  % keep old pairs and add from the new iterate   
            maskVg=1:numsy2;
            if numsy>0
                if Nflag==1  % quasi-Newton step was accepted
                    
                    %Vs=(-alpha/ns)*(Vg_old(maskVg)...
                    %    +VV(maskVV,maskVV)*p(1:numsy2));
                    
                    Vs = (V(:,maskV)'*(s/ns))./nV(maskV);
                    
                else
                    
                    %Vs=(-alpha/ns)*(Vg_old(maskVg)...
                    %   +VV(maskVV,maskVV(lindepV))*p(1:rankV));
                    
                    Vs = (V(:,maskV)'*(s/ns))./nV(maskV);
                    
                end
            end
            
            numsy   = numsy+1;
            numsy2  = numsy*2;
            maskV   = [1:numsy, m+1:m+numsy];
        else  % remove the oldest pair                        
            
            maskVc  = maskV([2:m m+2:m2]);            
            maskV   = maskV([ 2:m 1 m+2:m2 m+1]);
            maskVg  = [2:m m+2:m2];
            
            if Nflag==1
                
                %Vs=(-alpha/ns)*(Vg_old(maskVg)+VV(maskVg,:)*p(1:numsy2));
                Vs = (V(:,maskVc)'*(s/ns))./nV(maskVc);
            else
                
                %Vs=(-alpha/ns)*(Vg_old(maskVg)+VV(maskVg,lindepV)*p(1:rankV));                
                Vs = (V(:,maskVc)'*(s/ns))./nV(maskVc);                
            end
            
            E(1:m-1)                                = E(2:m);
            VV([1:m-1,m+1:m2-1],[1:m-1,m+1:m2-1])   = VV([2:m,m+2:m2],[2:m,m+2:m2]);
            T(1:m-1,1:m-1)                          = T(2:m,2:m);
            L(1:m-1,1:m-1)                          = L(2:m,2:m);           
        end         
        
        E(numsy)                    = sy/ny^2;            
        V(:,maskV(numsy))           = s;
        nV(maskV([numsy,numsy2]))   = [ns;ny];
        V(:,maskV(numsy2))          = y;        
        VV([numsy,m+numsy],numsy)   = [1; sy/ns/ny];    
        
        VA(maskV(numsy2),1:rankA)   = (Ag(1:rankA)-Ag_old(1:rankA))'./ny;
        VA(maskV(numsy),1:rankA)    = (A(maskA,:)*s)'./ns;
        
        if numsy>1            
            VV([1:numsy-1 m+1:m+numsy-1],numsy)=Vs;
        end
        
        VV([numsy,m+numsy],m+numsy)             = [sy/ns/ny;1];        
        Vg(1:numsy2)                            = (V(:,maskV)'*g)./nV(maskV);
        VV([1:numsy-1 m+1:m+numsy-1],m+numsy)   = (Vg([1:numsy-1,numsy+1:numsy2-1])...
            -Vg_old(maskVg))/ny;
        VV(numsy,[1:numsy-1 m+1:m+numsy-1])     = VV([1:numsy-1 m+1:m+numsy-1],numsy);
        VV(m+numsy,[1:numsy-1 m+1:m+numsy-1])   = VV([1:numsy-1 m+1:m+numsy-1],m+numsy);
        T(1:numsy,numsy)                        = VV(1:numsy,m+numsy); 
        L(numsy,1:numsy-1)                      = VV(numsy,m+1:m+numsy-1); 
        
        Vg_old(1:numsy2)                        = Vg(1:numsy2);            
        Ag_old(1:rankA)                         = Ag(1:rankA);
    else  % skip L-BFGS update but compute V'*g        
        Vg(1:numsy2)        =(V(:,maskV)'*g)./nV(maskV);
        Vg_old(1:numsy2)    = Vg(1:numsy2);   
        
        Ag(1:rankA)         = A(maskA,:)*g;
        Ag_old(1:rankA)     = Ag(1:rankA);
        
        numskip             = numskip + 1;
    end  % L-BFGS update
       
    %%--------------------------------- Equality constrained step ----------------------------%%
    % Computes the step from the formula (11) in the manuscript:
    % (Ki)_11 = 1/delta.I + [V A] M [V A]^T
        
    alpha   = 1/delta;
    gamma   = delta;
    
    Uc(numsy+1:numsy2,1:(rankA+1))  = linsolve(T(1:numsy,1:numsy),...
                                        [-VA(maskV(1:numsy),1:rankA) -Vg(1:numsy)],...
                                            opts1); % Triangular system
    
    Uc(1:numsy,1:(rankA+1))         = linsolve(T(1:numsy,1:numsy),...
                                        -((diag(gamma.*E(1:numsy)) +VV(m+1:m+numsy,m+1:m+numsy) )*...
                                        Uc(numsy+1:numsy2,1:(rankA+1))+...
                                        [VA(maskV(numsy+1:numsy2),1:rankA) Vg(numsy+1:numsy2)])... 
                                            ,opts2); % Transposed triangular system
                                        
    b1(1:rankA)         = (-gamma.*b(maskA) + (Ag(1:rankA)+ Uc(1:numsy2,1:rankA)'*Vg(1:numsy2) ));
    b2(1:rankA)         = linsolve( (AA(1:rankA,1:rankA) + Uc(1:numsy2,1:rankA)'*VA(maskV(1:numsy2),1:rankA)),...
                            b1(1:rankA),opts3); % Symmetric linear system
    
    gradp               = g-A(maskA,:)'*b2(1:rankA);
    
    Vgradp(1:numsy2)    = Vg(1:numsy2)-VA(maskV(1:numsy2),1:rankA)*b2(1:rankA);
    
    st                  = -alpha.*(V(:,maskV)*((Uc(1:numsy2,rankA+1)- Uc(1:numsy2,1:rankA)*b2(1:rankA))./nV(maskV)) + gradp);
    
    nst                 = norm(st);    
    
    af=max(1,abs(f));
    if nst<= (trrad+trradb)  % quasi-Newton step is inside TR
        
        % Compute the trial point and evaluate function    
        xt      = x+st;
        [ft,gt] = fun(xt);
        numf    = numf+1;
        numg    = numg+1;
        
%         [A,b]   = const(xt);
%         numc    = numc+1;
%         numcg   = numcg+1;
        
        % Compare actual (ared) and predicted (pred) reductions 
        ared = ft-f;
        if abs(ared)/af<ftol
              % relative actual reduction is too small

            ratio = 1;
            
        elseif ared<0
            %pred=-0.5*alpha*(ng^2+p(1:numsy2)'*Vg(1:numsy2));
            %ratio=ared/pred;                
            
            % This version accepts any step that reduces the objective function value            
            ratio = 1;
        else
            ratio = 0;
        end
        
        if ratio>tau0                   
            Nflag = 1;  % quasi-Newton step is accepted    
        else
            Nflag = 0;  % quasi-Newton step is rejected by ratio
        end
        
        if (ratio<tau1)
            trrad=min(c1*nst,c2*trrad);
            %trrad = max(c1*nst,c2*trrad);
        end   
        
        if trrad<trrad_min  % TR radius is too small, terminate
            ex      = -1;
            tcpu    = toc(t0);
            outinfo = LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,numskip,numgp0);
            return
        end
                    
    else  % quasi-Newton step is rejected by TR
        Nflag=0;  
    end  % checking quasi-Newton step
    
    if ~Nflag
        %% ----- Quasi-Newton step is rejected, Newton's method ----- %%
        
        tract   = tract+1;
        trit    = 0;
        ratio   = 0;
            
        while (ratio<=tau0)
            
            itroot  = 0;
            phi     = 1/nst - 1/trrad;
            
            % Newton's method
            sigma   = 0;
                        
            while abs(phi) > epsroot && itroot < maxitroot
                
                gam     = gamma+sigma;
                del     = 1/gam;
                alp     = gamma/gam;
                a1      = 1/alp-1;
                
                
                N12(1:numsy,1:numsy)    = a1.*VV(1:numsy,m+1:m+numsy)-T(1:numsy,1:numsy)./alp;
                Uc(1:numsy2,1:rankA+1)  = linsolve([a1.*VV(1:numsy,1:numsy),N12(1:numsy,1:numsy); N12(1:numsy,1:numsy)',...
                                            -(diag(gam.*E(1:numsy))+VV(m+1:m+numsy,m+1:m+numsy))],...
                                            [VA(maskV(1:numsy2),1:rankA),Vg(1:numsy2)],opts3);
                if itroot > 0

                    b1(1:rankA)         = (-gam.*b(maskA) + (Ag(1:rankA)+ Uc(1:numsy2,1:rankA)'*Vg(1:numsy2) ));
                    b2(1:rankA)         = linsolve( (AA(1:rankA,1:rankA) + Uc(1:numsy2,1:rankA)'*VA(maskV(1:numsy2),1:rankA)), b1(1:rankA), opts3); % Symmetric linear system

                    gradp               = g-A(maskA,:)'*b2(1:rankA);

                    Vgradp(1:numsy2)    = Vg(1:numsy2)-VA(maskV(1:numsy2),1:rankA)*b2(1:rankA);

                    st                  = -del.*(V(:,maskV)*((Uc(1:numsy2,rankA+1)- Uc(1:numsy2,1:rankA)*b2(1:rankA))./nV(maskV)) + gradp);
                        
                end
                
                Vs(1:numsy2)            = (V(:,maskV)'*st)./nV(maskV);
                As(1:rankA)             = A(maskA,:)*st;
                Uc(1:numsy2,rankA+1)    = linsolve([a1.*VV(1:numsy,1:numsy),N12(1:numsy,1:numsy); N12(1:numsy,1:numsy)',...
                                            -(diag(gam.*E(1:numsy))+VV(m+1:m+numsy,m+1:m+numsy))],...
                                            Vs(1:numsy2),opts3);
                                        
                b1(1:rankA)         = (-gam.*b(maskA) + (As(1:rankA)+ Uc(1:numsy2, 1:rankA)'*Vs(1:numsy2) ));
                b2(1:rankA)         = linsolve( (AA(1:rankA,1:rankA) + Uc(1:numsy2,1:rankA)'*VA(maskV(1:numsy2),1:rankA)), b1(1:rankA), opts3); % Symmetric linear system
    
                stp                 = -del.*(V(:,maskV)*((Uc(1:numsy2,rankA+1)- Uc(1:numsy2,1:rankA)*b2(1:rankA))./nV(maskV)) + st-A(maskA,:)'*b2(1:rankA));
                
                nst                 = norm(st);
                phi                 = 1/nst - 1/trrad;
                phip                =-(stp'*st)/nst^3;
                

                sigma   = sigma - phi/phip;
   
                itroot  = itroot+1;
                            
            end

            xt      =x+st;
            ft      = fun(xt);
            numf    = numf+1;
            
            % Check for function reduction
            ared=ft-f;                        
            if abs(ared)/af<ftol
                ratio=1;
            elseif ared<0               
                %pred=vpar(1:rankV)'*(gpar(1:rankV)+0.5*lam(1:rankV).*vpar(1:rankV))...
                %    +(alpha^2*gp/2-alpha)*ngperp^2;                
                %ratio=ared/pred;
                
                % This version accepts any step that reduces the function
                % value.
                ratio = 1;
            else
                ratio = 0;
            end              
                        
            if (ratio<tau1)
                trrad=min(c1*nst,c2*trrad);
                %trrad = max(c1*nst,c2*trrad);
            end              
                        
            if trrad<trrad_min  % TR radius is too small, terminate
                ex      = -1;
                tcpu    = toc(t0);
                outinfo = LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,numskip,numgp0);
                return
            end 
            
            trit = trit+1;
        end  % solving the TR subproblem
        
        % Update counter if the initial trial step is rejected and TR
        % subproblem is to be solved again for a reduced TR radius
        if trit>1
            trrej = trrej+1;
        end
    end    
    %% Trial step is accepted, update TR radius and gradient
    if (ratio>tau2) && (nst>=c3*trrad)   
        trrad = c4*trrad;   
    end                
    s   = st;
    ns  = norm(s);
    x   = xt;
    f   = ft;    
    g_old = g;
    if Nflag==1
        g = gt;
    else
        [~, g]      = fun(x);
        %numf        = numf+1;
        numg        = numg+1;
    end   
        
    ng  = norm(g);
    it  = it+1; 
    
    [~,b]           = const(x);    
    numc            = numc+1;
    nb              = norm(b(maskA));
    
    Ag(1:rankA)     = A(maskA,:)*g;        
    b2(1:rankA)     = linsolve(RA(1:rankA,1:rankA),...
                                    linsolve(RA(1:rankA,1:rankA),Ag(1:rankA),opts2),...
                                    opts1);
    
    % Computing \| Proj(g) \|^2
    ngp             = sqrt(abs(ng^2 - Ag(1:rankA)'*b2(1:rankA)));
    
    conv            = ngp/max(1,norm(x));
    
    % Display information about the last iteration
    if dflag==1 
        fprintf('%d\t %.3e\t %.3e\t %.3e\t %.3e\t %.3e\n',...
            it,f,nb,ngp,ns,trrad);                
    end    
    
    % Check the main stopping criteria
    if (conv>gtol) || (nb > btol)  % ng>gtol*max(1,norm(x))   
        y   = g-g_old;
        ny  = norm(y);
        sy  = s'*y;
    else
        mflag = 0;
    end
    
    % Check if the algorithm exceeded maximum number of iterations
    if it > maxit           
        ex      = -2;
        tcpu    = toc(t0);
        outinfo = LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,numskip,numgp0);
        return;
    end   
    
    if storedat ~= 0
         trrads(it)      = trrad; % store Delta_k,
         funs(it)        = f; % f_k,
         converr(it)     = conv;
    end
        
end  % main loop

ex=1;
tcpu=toc(t0);  
outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,numskip,numgp0);
