% RSQP Optimization Toolbox
% Version 1.0 (R14SP3) 10-Jul-2006 
%
% Nonlinear minimization of functions.
%   rsqp         - Multidimensional equality constrained nonlinear minimization.
%
% Auxiliary function.
%   nlconst_rsqp - Find the constrained minimum of a function. Called by rsqp.
%   qusub2       - Find minimum of quadratic programming problems.
%
% Controlling defaults and options.
%   optimset2    - Create or alter optimization OPTIONS structure. 
%   optimget2    - Get optimization parameters from OPTIONS structure. 
%
% Interface.
%   proxy        - Interface of RSQP Toolbox and Optimization Toolbox.
%
%  Options setting examples from User's Guide.
%   ex_bt12          - nonlinear programming problem bt12 from CUTE test set
%   ex_bt12_fun      - nonlinear objective with gradient of problem bt12
%   ex_bt12_con      - nonlinear constraints with gradients of problem bt12
%   ex_hs006         - nonlinear programming problem hs006 from CUTE test set
%   ex_hs006_fun     - nonlinear objective with gradient of problem hs006
%   ex_hs006_con     - nonlinear constraints with gradients of problem hs006
%   ex_lsnnodoc      - nonlinear programming problem lsnnodoc from CUTE test set
%   ex_lsnnodoc_fun  - nonlinear objective with gradient of problem lsnnodoc
%
% Large-scale examples from User's Guide
%   ex_halffree      - nonlinear programming problem halffree from Biegler et al.(2000)
%   ex_halffree_fun  - nonlinear objective with gradient of problem halffree
%   ex_halffree_con  - nonlinear constraints with gradients of problem halffree
%   ex_lalee         - nonlinear programming problem lalee from Lalee (1992)
%   ex_lalee_fun     - nonlinear objective with gradient of problem lalee
%   ex_lalee_con     - nonlinear constraints with gradients of problem lalee
%   ex_onefree       - nonlinear programming problem onefree from Biegler et al.(2000)
%   ex_onefree_fun   - nonlinear objective with gradient of problem onefree
%   ex_onefree_con   - nonlinear objective with gradient of problem onefree
%

%   Copyright 2006 Institute of Industrial Control, Zhejiang University.
%   $Date: 2006/07/10 13:03:26 $
