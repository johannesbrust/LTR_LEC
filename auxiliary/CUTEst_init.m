function CUTEst_init()

if ispc
  disp('Running on Windows');
  disp('No installation available!');
end

if isunix,
  disp('Running on Unix');
  
  % adapt the next two lines for your CUTEr installation folder
  % J.B, 10/13 adapted for CUTEst installation
  %addpath('/sw/opt/CUTEr/130204-r150/cuter2/CUTEr.large.pc.lnx.gfo/bin');
  addpath('~/Documents/CUTEST/cutest/src/matlab');
end