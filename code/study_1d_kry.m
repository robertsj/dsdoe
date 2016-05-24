% 1-D, 1-G study for 18086 project
%  -- S8, sigA = sigS 0.5, S = 1
%  -- 10 cm slab
%  -- mesh count ranges = [ 10 20 50 100 250 500 1000 2000 3000 4000 5000]
%  uses PI, gauss elimination, gmres+ilu(0), gmres+ilu(1e-4), 
%    bicg-stab+ilu(0), andbicg-stab+ilu(1e-4)
%  tracks computation time for each method; breaks krylov solver time into
%    K matrix, LU, and solver time; for elim, just K matrix and solver

clear 
format short
mesh=1;
% set of mesh counts
t_pi     = zeros(length(mesh),1);
t_matx   = zeros(length(mesh),1);
t_elim   = zeros(length(mesh),1);
t_ilu0   = zeros(length(mesh),1);
t_iluT   = zeros(length(mesh),1);
t_gmres0 = zeros(length(mesh),1);
t_gmresT = zeros(length(mesh),1); 
t_bicg0  = zeros(length(mesh),1);
t_bicgT  = zeros(length(mesh),1); 

for mm = 10000:10000:50000

data    = [ 1.0 0.5 0.5
            0.4 0.3 0.1
            2.0 0.5 1.5 ];
xcm     = [ 0   10    20   30 ];
xfm    = [  1   1    1     ]*mm;
numg    = 1;
mt      = [  1 2 3 ];
src     = [ 1 0  2 ];

input   =   struct(   ...
    'numg',         numg, ...     % number of groups
    'numm',            2, ...     % number of materials
    'xcm',           xcm, ...     % slab bounds
    'xfm',           xfm, ...     % number of fine meshes
    'mt',             mt, ...     % slab material ids
    'data',         data, ...     % mat comp's
    'src',           src, ...     % volume source
    'ord',             8, ...     % number of ordinates
    'maxit',        1000, ...     % max iterations
    'maxerr',       1e-6, ...     % max pointwise phi error
    'adj',             0, ...     % adjoint flag
    'bcL',             0, ...
    'bcR',             0 ...
    );

%-----MATRIX FORMULATION

tic
[KK,Q] = sn_one_d_2g_matrix_vec(input);
t_matx = toc;

for i = 1:length(mesh)

    %---DIRECT SOLVERS

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
disp([' matrix time     : ',num2str(t_matx)])
disp([' time elimination: ',num2str(t_elim)])
disp([' time ilu(0)     : ',num2str(t_ilu0)])
disp([' time ilu(t)     : ',num2str(t_iluT)])
disp([' time gmres(0)   : ',num2str(t_gmres0), ' tot = ',num2str(t_gmres0+t_ilu0)])
disp([' time gmres(T)   : ',num2str(t_gmresT), ' tot = ',num2str(t_gmresT+t_iluT)])
disp([' time bicg(0)    : ',num2str(t_bicg0), ' tot = ',num2str(t_bicg0+t_ilu0)])
disp([' time bicg(T)    : ',num2str(t_bicgT), ' tot = ',num2str(t_bicgT+t_iluT)])
%time_plot(mesh,t_pi,t_matx,t_elim,t_ilu0,t_iluT,t_gmres0,t_gmresT,t_bicg0 ,t_bicgT)
%title('1-D Mesh Study')
end