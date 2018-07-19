function [x,f,outinfo] = LTRSC_LEC(fun,const,x0,params)
% LTRSC_LEC - Limited-Memory Trust-Region Shape-Changing Norm for
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
% Here || s ||_{P,infty} is the shape-changing norm, Delta_k 
% is the "trust-region" radius,  and B = B_k is the L-BFGS
% compact representation: 
%
%       B = delta.I + V \hat{Xi} V', and B^{-1} = H = gamma.I - V \hat{M} V'. 
%
% The notation corresponding to the article is V (code) = \hat{\Psi} (article) and 
% delta (code) = gamma (article), \hat{M} (code) = -\hat{M} (article).
% The initial parameters are defined in a struct params. 
% The output information is contained in a struct outinfo. 
%
% params contains the following fields {default values}:
%   m       - maximum number of stored vector pairs (s,y) {5}
%   [gtol    - tolerance on L2-norm of gradient ||Proj(g)||<gtol*max(1,||x||)
%   {1e-5}]
%   ftol    - tolerance for abs(f_k-f_{k+1})
%   maxit   - maximum number of iterartions {100000}
%   ranktol - tolerance threshold for establishing rank of V {1e-7}
%   dflag   - display parameter {1}:
%               1 - display every iteration;
%               0 - no display.
%   storedat- flag whether or not to store data.
%               1 - store:
%                       - Delta_k,                       
%                       - f_k
%                       - Convergence criterion,
%   trradb  - trrad_k = trradb + trrad (for linearly constrained
%   problems)
%   btol    - ||h_k - b_h||<btol
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

J.B 20/01/17 
Modification to allow for different choices of gamma_perp in the
trust-region subproblem. This is only applicable if the quasi-Newton step
is rejected. The L-BFGS matrix has form

    B_k = B_0 + V M V',

             | gamma.I          0     |
    B_0 =  P |   0        gamma_perp.I| P'

J.B 08/23/17
Apply shape-changing norm to 

	||Wg||_{P,inf} 	instead of ||g||_{P,inf},

	where

	W = I - A'(A B^{-1} A')^{-1} A B^{-1} (oblique projection)

J.B 08/24/17
Modification of stopping criterion. Instead of the true solution, seek a
point for which || Proj(g) || = 0.

J.B 10/16/17
Compact representation of the solution matrix

    B^{-1} W = delta.I - Psih Mbb Psih',

where Psih = [A' V] and the matrix Mbb (2m+mm x 2m+mm) is composed as

    Mbb = [ delta.Mb  Mb*Uct; (Mb*Uct)' -N + Uct*V*A'],

where Uct = A*V*M. 

This approach computes the eigendecomposition of B^{-1} W and applies the
shape-changing norm to the vector B^{-1} W g.

J.B., 11/10/17
Computing eigenvalues of the compact representation using a deflation
technique.

    B^{-1} W = delta.I - Psih Mbb Psih'
             

            = [P_||1  P_||2 P_|_] *
                |0      0           0|
                |0  delta.I - D     0| *
                |0      0     delta.I|   
              [P_||1  P_||2 P_|_]'.

J.B., 06/27/18
Explicit treatment of linear matrix A. Sparse capabilities.

J.B., 07/02/18 
Efficient updating of AV only. Inclusion of efficient convergence test, by
using a LDL' factorization of AA'. This allows rank deficient matrices A,
too.

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
    ranktol = 1e-10; % 1e-7
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
t0 = tic;     

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

gp_old  = 0; % Gamma for dense B0.

numc    = 0; % number of constraint evaluations
%numcg   = 0; % number of constraint gradient evaluations 

% Constraints

[A,b]       = const(x0);
numc        = numc + 1; % Internal counter constraint evaluations.
%numcg       = numcg + 1; % Internal counter constraint gradient evaluations.
mm          = length(b);
maskA       = 1:mm;
rankA       = length(maskA);

AA          = zeros(rankA,rankA); % AA'

RA          = zeros(mm,mm); % RA'*RA = AA', Cholesky or LDL' factor of AA'

VA          = zeros(m2,mm); % V'A'

Uc          = zeros(m2,mm+1); % [gamma.* N] [ V'A' V'g ]
b1          = zeros(mm,1); % b+...
b2          = zeros(mm,1); % (AA'+...)^{-1} b2
Ag          = zeros(mm,1); % A*g
Ag_old      = zeros(mm,1); % A*g previous iteration
Vgradp      = zeros(m2,1); % V'*W*g

% Projected compact representation
ml      = m2+mm;

alpha   = 0; %#ok<NASGU>
p       = zeros(ml,1);
 
gpar    = zeros(ml,1);  % gpar=Ppar*g, where Ppar=inv(R)*V*U
agpar   = zeros(ml,1);  % agpar=abs(gpar)
vpar    = zeros(ml,1);  % vpar=Ppar*s

Ti      = zeros(m,1); % inv(T)
N       = zeros(ml,ml);  % N=[inv(T)'(E+delta.Y'Y)inv(T) -delta.inv(T)'; -delta.inv(T) 0]
Mb      = zeros(mm,mm); % Mb = delta.(AA' + Uct*V'*A')^{-1}

R       = zeros(ml,ml);  % upper-triangular matrix in QR decomposition of Psih = [A' V]
U       = zeros(ml,ml);  % orthogonal matrix, eigenvectors of R*W*R'=U*D*U'
D       = zeros(ml,ml);  % diagonal matrix, eigenvalues of R*W*R'=U*D*U'
lam     = zeros(ml,1);  % vector, lam = diag(delta.I-D);

MbUct   = zeros(mm,m2); % Mb*Uct

AAV     = zeros(mm,m2+1); % (AA')^{-1} [VA' Ag]
p2      = zeros(mm,1); % (AA')^{-1} p

% Precomputed Identities, and zero matrix
Im      = eye(m); 
Imm     = eye(mm); 
zm      = zeros(m,m); 

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
idxnzr  = find(sA > ranktol);
dsA     = diag(sA);
% Upper triangular matrix, and rank determination of A.
[Ldla, Ddla, pdla]  = ldl(dsA(idxnzr,(idxnzr))\...
                        ((A((idxnzr),:)*A((idxnzr),:)')/...
                        dsA((idxnzr),(idxnzr))),'vector');
                    
dDdla               = diag(Ddla);

% Compute linearly independent columns.
maskldla            = find(dDdla>1e-10); % ^2
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
ngp                     = sqrt(abs(ng*ng-Ag(1:rankA)'*b2(1:rankA)));
d                       = d_./ngp;

conv                    = ngp/max(1,norm(x0));

if dflag==1
    fprintf('\n**********************\nRunning LTRSC-LEC\n');
    fprintf('it\t obj\t\t norm(b)\t norm(Proj(g)) \t norm(dx)\t trrad\n');
    fprintf('%d\t %.3e\t %.3e\t %.3e\t---\t\t ---\n',0,f0,nb,ngp);
end

if ((conv<gtol) && (nb < btol)) % ng<max(1,norm(x0))*gtol,
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

ngp             = sqrt(abs(ng*ng-Ag(1:rankA)'*b2(1:rankA)));
conv            = ngp/max(1,norm(x));

[~,b]           = const(x);
numc            = numc+1;

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
            
            outinfo = LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,numskip,numgp0);
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
                  %Vs = (V(:,maskV(lindepV))'*(s/ns))./nV(maskV(lindepV));
                  
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
            VV([1:numsy-1 m+1:m+numsy-1],numsy) = Vs;
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
        
        ared=ft-f;
        if abs(ared)/af<ftol
              % relative actual reduction is too small

            ratio =1;
            
        elseif ared<0
            %pred=-0.5*alpha*(ng^2+p(1:numsy2)'*Vg(1:numsy2));
            %ratio=ared/pred;            
            
            %This version accepts any step that reduces the objective function value
            ratio =1;
        else
            ratio=0;
        end
        
        if ratio>tau0                   
            Nflag=1;  % quasi-Newton step is accepted
            nst_inf_2=nst;      
        else
            Nflag=0;  % quasi-Newton step is rejected by ratio
        end
        
        if (ratio<tau1)
            trrad=min(c1*nst,c2*trrad);
            %trrad=max(c1*nst,c2*trrad);
        end   
        
        if trrad<trrad_min  % TR radius is too small, terminate
            ex=-1;
            tcpu=toc(t0);
            outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,numskip,numgp0);
            return
        end
                    
    else  % quasi-Newton step is rejected by TR
        Nflag=0;  
    end  % checking quasi-Newton step
    
    if ~Nflag
        %%------------------ Eigenvalues of Bi W ---------------------%%
        %% 10/17/17 Eigendecomposition of compact representation of matrix         
        % Version to compute eigenvectors in P_||2 only.
        
        tract   = tract+1;        
        idelta  = 1/delta;
        
        AAV(1:rankA,1:numsy2 +1) = linsolve(RA(1:rankA,1:rankA),...
                                    linsolve(RA(1:rankA,1:rankA),...
                                        [VA(maskV(1:numsy2),1:rankA)' Ag(1:rankA)],opts2),...
                                            opts1);       
        %AAV(1:rankA,1:numsy2 +1) = linsolve(AA(1:rankA,1:rankA), [VA(maskV(1:numsy2),1:rankA)' Ag(1:rankA)], opts3);
        
        % Compute LDL decomposition of VV(pdl,pdl)=Ldl*Ddl*Ldl'
        % and R=sqrt(Ddl)*Ldl' in the QR decomposition of Psi = [A' V].
        [Ldl, Ddl, pdl] = ldl(VV([1:numsy,m+1:m+numsy],[1:numsy,m+1:m+numsy])-...
                                VA(maskV(1:numsy2),1:rankA)*AAV(1:rankA,1:numsy2),'vector');  
        dDdl            = diag(Ddl);
        
        % Use only safely linearly independent columns of Psi to
        % compute the trial step
        maskldl         = find(dDdl>ranktol); % ^2 
        rankV           = length(maskldl);  % column rank of Psi
        lindepV         = pdl(maskldl);  % index of safely linearly independent columns of Psi
        
        % Compute inverse permutation of pdl
        ipdl                = 1:numsy2;          
        ipdl(pdl)           = ipdl;
                 
        R(1:rankV,1:numsy2) = diag(sqrt(dDdl(maskldl)))*Ldl(:,maskldl)';
        
        % Compute the lower-right hand block of the middle matrix Mbb in HW = gamma*I-Psi*M*Psi'
        Ti(1:numsy,1:numsy)     = linsolve(T(1:numsy,1:numsy),Im(1:numsy,1:numsy),opts1);
        N(1:numsy2,1:numsy2)    = [Ti(1:numsy,1:numsy)'*(diag(E(1:numsy))+...
            idelta.*VV(m+1:m+numsy,m+1:m+numsy))*Ti(1:numsy,1:numsy) -idelta.*Ti(1:numsy,1:numsy)';...
            -idelta.*Ti(1:numsy,1:numsy) zm(1:numsy,1:numsy)];
        
        Mb(1:rankA,1:rankA)     = idelta.*(linsolve((AA(1:rankA,1:rankA)+Uc(1:numsy2,1:rankA)'*VA(maskV(1:numsy2),1:rankA)),...
                                Imm(1:rankA,1:rankA),opts3));
                            
        MbUct(1:rankA,1:numsy2) = Mb(1:rankA,1:rankA)*Uc(1:numsy2,1:rankA)';
        
        % Middle matrix with Triangular matrices around M1 = R M R'
        M1                      = R(1:rankV,ipdl)*...
                                (Uc(1:numsy2,1:rankA)*MbUct(1:rankA,1:numsy2)-N(1:numsy2,1:numsy2))*R(1:rankV,ipdl)';
        
        [U(1:rankV,1:rankV),D(1:rankV,1:rankV)] = eig(M1);        
        lam(1:rankV)                            = idelta-diag(D(1:rankV,1:rankV));
                 
        gpar(1:rankV)                           = U(1:rankV,1:rankV)'*linsolve(R(1:rankV,maskldl),(Vg(lindepV) - VA(maskV(lindepV),1:rankA)*...
                                                    AAV(1:rankA,numsy2+1)),opts2);
                                                
        agpar(1:rankV)  = abs(gpar(1:rankV));
        ngpar           = norm(gpar(1:rankV));
        ngperp          = sqrt(max(0,(ng-ngpar)*(ng+ngpar)));    
        if (ngperp*ngperp)/(ng*ng) < 1e-13; numgp0 = numgp0+1; end;
        
        vpar(1:rankV)   = -gpar(1:rankV).*lam(1:rankV); 
        trit=0;
        
        %% Solving the TR subproblem in the new variables
        % Corresponds to eqs. (28,29) in the article
        ratio=0;
               
        % Takes the running maximum of in the space of P_{\perp} of 
        % the eigendecomposition 
        % (Ki) = [P P_{\perp}] diag(lam, gp) [P P_{\perp}]^T.
        % The matrices P, P_{\perp} are not explicitly computed in the
        % code.
        if gp_old < delta; gp_old = delta; end  
        gp = gp_old;

        while (ratio<=tau0)
                        
            % Compute alpha and vpar
            alpha   = 1/gp;
            if 1e-14 < abs(ngperp)
                alpha = min(trrad/ngperp, 1/gp);
            end
            %if (1/gp < (trrad)/ngperp); alpha=1/gp; else; alpha = (trrad)/ngperp; end
            %alpha=min(1/gp,(trrad)/ngperp);
            
            ns_perp     = min(trrad,ngperp/gp);  % L2 norm of s_perp
            ind         = find(lam(1:rankV).*agpar(1:rankV)>(trrad));
            vpar(ind)   = -trrad*sign(gpar(ind));
            nspar       = norm(vpar(1:rankV),inf);  % L2 norm of spar
            nst_inf_2   = max(nspar,ns_perp);  % new (Linf-L2) norm of spar
            
            p(1:rankV)  = linsolve(R(1:rankV,maskldl),...
                    U(1:rankV,1:rankV)*(vpar(1:rankV)+alpha*gpar(1:rankV)),opts1);
            
            p2(1:rankA) = AAV(1:rankA,lindepV)*p(1:rankV);
            
            % Forming trust-region step
            st          = -alpha*g+ A(maskA,:)'*(alpha*AAV(1:rankA,numsy2+1)-p2(1:rankA)) ...
                            +V(:,maskV(lindepV))*(p(1:rankV)./nV(maskV(lindepV)));            
            %st=-alpha*(g+V(:,maskV(lindepPsi))*(p(mm+1:numml)./nV(maskV(lindepPsi)))+A(lindepA,:)'*p(1:mm));

            xt      = x+st;
            [ft]    = fun(xt);
            numf    = numf+1;
            
            % Compare actual (ared) and predicted (pred) reductions 
            ared = ft-f;                        
            if abs(ared)/af<ftol
                ratio = 1;
            elseif ared<0               
                %pred=vpar(1:rankV)'*(gpar(1:rankV)+0.5*lam(1:rankV).*vpar(1:rankV))...
                %    +(alpha^2*gp/2-alpha)*ngperp^2;                
                %ratio=ared/pred;
                
                % This version does not compute the predicted reduction
                ratio = 1;
            else
                ratio = 0;
            end              
                        
            if (ratio<tau1)
                trrad = min(c1*nst_inf_2,c2*trrad);
                %trrad=max(c1*nst,c2*trrad);
            end              
                        
            if trrad<trrad_min  % TR radius is too small, terminate
                ex = -1;
                tcpu = toc(t0);
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
    if (ratio>tau2) && (nst_inf_2>=c3*trrad)   
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
