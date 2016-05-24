% 1-D, 1-G study for 18086 project
%  -- S8, sigA = sigS 0.5, S = 1
%  -- 10 cm slab
%  -- mesh count ranges = [ 10 20 50 100 250 500 1000 2000 3000 4000 5000]
%  uses PI, gauss elimination, gmres+ilu(0), gmres+ilu(1e-4), 
%    bicg-stab+ilu(0), andbicg-stab+ilu(1e-4)
%  tracks computation time for each method; breaks krylov solver time into
%    K matrix, LU, and solver time; for elim, just K matrix and solver

clear, clc 
format short

% set of mesh counts
sn       = [ 2 4 8 16 32];
t_pi     = zeros(length(sn),1);
t_matx   = zeros(length(sn),1);
t_elim   = zeros(length(sn),1);
t_ilu0   = zeros(length(sn),1);
t_iluT   = zeros(length(sn),1);
t_gmres0 = zeros(length(sn),1);
t_gmresT = zeros(length(sn),1); 
t_bicg0  = zeros(length(sn),1);
t_bicgT  = zeros(length(sn),1); 

src     = [ 1 ];
data    = [ 1.0 0.5 0.5];
xcm     = [ 0   10 ];
xfm     = [ 2000 ];
numg    = 1;
mt      = [  1  ];

for i = 1:length(sn)
    disp(['sn: ',num2str(sn(i))])
    input   =   struct(   ...
        'numg',         numg, ...     % number of groups
        'numm',            1, ...     % number of materials
        'xcm',           xcm, ...     % slab bounds
        'xfm',           xfm, ...     % number of fine meshes
        'mt',             mt, ...     % slab material ids
        'data',         data, ...     % mat comp's
        'src',           src, ...     % volume source
        'ord',         sn(i), ...     % number of ordinates
        'maxit',        1000, ...     % max iterations
        'maxerr',       1e-6, ...     % max pointwise phi error
        'adj',             0, ...     % adjoint flag
        'bcL',             0, ...
        'bcR',             0 ...
        );

    %-BEGIN SOLVERS HERE
    
    %---POWER ITERATION
    tic
    [phi,psi]   = sn_one_d(input);
    t_pi(i) = toc;
    clear phi psi
    pause(5)

    %---DIRECT SOLVERS
    
    %-----MATRIX FORMULATION

    tic
    [KK,Q] = sn_one_d_2g_matrix_vec(input);
    t_matx(i) = toc;
    
    
    %-----GAUSSIAN ELIMINATION
    tic
    p = KK\Q;
    t_elim(i) = toc;
    clear p
    pause(5)

    %-----ILU(0)
    setup.type='nofill';
    tic
    [L,U] = ilu(KK,setup);
    t_ilu0(i) = toc;
    
    %-------GMRES+ILU(0)
    tic
    p=gmres(KK,Q,4,input.maxerr,100,L,U);
    t_gmres0(i) = toc;
    clear p
    pause(5)
    
    %-------BICG-STAB+ILU(0)
    tic
    p=bicgstab(KK,Q,input.maxerr,100,L,U);
    t_bicg0(i) = toc;
    clear p L U
    pause(5)
    
    %-----ILU(1e-4)
    setup.type = 'ilutp';
    setup.droptol = 0.0001;
    tic;
    [L,U] = ilu(KK,setup);
    t_iluT(i) = toc;
    
    %-------GMRES+ILU(1e-4)
    tic
    p=gmres(KK,Q,4,input.maxerr,100,L,U);
    t_gmresT(i) = toc;
    clear p
    pause(5)
    
    %-------BICG-STAB+ILU(1e-4)
    tic
    p=bicgstab(KK,Q,input.maxerr,100,L,U);
    t_bicgT(i) = toc;
    clear p L U
    pause(5)
    
end

time_plot(sn,t_pi,t_matx,t_elim,t_ilu0,t_iluT,t_gmres0,t_gmresT,t_bicg0 ,t_bicgT)
title('1-D S_N Order Study')
figure(1),xlabel('S_N order')