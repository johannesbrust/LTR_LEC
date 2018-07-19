% Extending function for parallel computation
function [ex,numf,numg,numit,tcpu,tract,numrst,numskip,var1,var2,var3,var4]=runAlgorithm_const_quad(alg,...
        x0,params,numRuns,varargin)
tcpu=zeros(numRuns,1);

if nargin == 4
    for i=1:numRuns
        [~,~,outinfo]=alg(@cutest_fun,@const_quad_2,x0,params);
        %[~,~,outinfo]=alg(@cuter_fun,x0,params);
        tcpu(i)=outinfo.tcpu;
    end
elseif nargin == 5
    for i=1:numRuns
        [~,~,outinfo]=alg(@cutest_fun,@cutest_const,x0,params);
        %[~,~,outinfo]=alg(varargin{1},x0,params);
        %[~,~,outinfo]=alg(@cuter_fun,x0,params);
        tcpu(i)=outinfo.tcpu;
    end
end
numf=outinfo.numf;
numg=outinfo.numg;
ex=outinfo.ex;
numit=outinfo.numit;   
tract=outinfo.tract;
if nargout==7
    numrst=outinfo.numrst;
elseif nargout==8
    numrst=outinfo.numrst;
    numskip=outinfo.numskip;
elseif nargout == 9
    numrst  =outinfo.numrst;
    numskip =outinfo.numskip;
    var1    = outinfo.varout2;
elseif nargout == 10
    numrst  =outinfo.numrst;
    numskip =outinfo.numskip;
    var1    = outinfo.varout2;
    var2    = outinfo.varout3;
elseif nargout == 11
    numrst  =outinfo.numrst;
    numskip =outinfo.numskip;
    var1    = outinfo.varout2;
    var2    = outinfo.varout3;
    var3    = outinfo.varout4;
elseif nargout == 12 
    numrst  =outinfo.numrst;
    numskip =outinfo.numskip;
    var1    = outinfo.varout2;
    var2    = outinfo.varout3;
    var3    = outinfo.varout4;
    var4    = outinfo.varout5;
end
