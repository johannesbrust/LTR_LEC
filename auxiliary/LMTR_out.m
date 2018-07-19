function outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst,varargin)
%   LMTR_out returns output information from LMTR:
%           ex      - exitflag:
%                            1, if  norm of gradient is too small
%                           -1, if  TR radius is too small, trrad<1e-15
%                           -2, if  the max number of iterations is exceeded
%           numit   - number of succesful TR iterations
%           numf    - number of function evaluations
%           numg    - number of gradient evaluations
%           tcpu    - CPU time of algorithm execution
%           tract   - number of iterations when TR is active
%           trrej   - number of iterations when initial trial step is rejected
%           params  - input paramaters
%           numrst  - number of restarts, optional
%           varargin- optional
%               varargin{1} - number of update skips, when s'y <= 0.
outinfo=struct;
outinfo.numf=numf;
outinfo.numg=numg;
outinfo.tcpu=tcpu;
outinfo.ex=ex;
outinfo.numit=it;
outinfo.tract=tract;
outinfo.trrej=trrej;
outinfo.params=params;
if nargin==9
    outinfo.numrst=numrst;
elseif nargin ==10
    outinfo.numrst=numrst;
    outinfo.numskip=varargin{1};
elseif nargin ==11
    outinfo.numrst=numrst;
    outinfo.numskip=varargin{1};
    outinfo.varout2=varargin{2};
elseif nargin ==12
    outinfo.numrst=numrst;
    outinfo.numskip=varargin{1};
    outinfo.varout2=varargin{2};
    outinfo.varout3=varargin{3};
elseif nargin ==13
    outinfo.numrst=numrst;
    outinfo.numskip=varargin{1};
    outinfo.varout2=varargin{2};
    outinfo.varout3=varargin{3};
    outinfo.varout4=varargin{4};
elseif nargin ==14
    outinfo.numrst=numrst;
    outinfo.numskip=varargin{1};
    outinfo.varout2=varargin{2};
    outinfo.varout3=varargin{3};
    outinfo.varout4=varargin{4};
    outinfo.varout5=varargin{5};
else
end