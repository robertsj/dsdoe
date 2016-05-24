% This input tests the 2-D matrix formulation and gives a matrix
%  like that of  Figure 6(b)
% j. roberts

% input contains
%     numg   = number of energy groups
%     numm   = number of materials
%     xcm    = coarse divisions
%     xfm    = fine mesh interval for coarse divisions
%     mt     = material assignment for each coarse division
%     data   = cross-sections in the form
%              (mat1/g1) sigTOT   sigA  sigSg1->g1  sigSg2->g1 ...
%                 ...g2) sigTOT   sigA  sigSg1->g2  ....
%              (mat2/g1) ...
%     src    = volumetric (isotropic) source by coarse mesh and energy
%              (for adj, src is the cross-section for the source region)
%     ord    = number of ordinates (2,4,8, or 12)
%     maxit  = maximum iterations
%     maxerr = maximum relative pointwise error in phi
clear

format short

xcm    = [ 0  1  2];
xfm    = [  1   1]*2;
ycm    = [ 0   1  2];
yfm    = [  1   1]*1;
src(1,:,:)    = [   1 1; 1 1 ]' ;
src(2,:,:)    = [   0 0; 0 0  ] ;                  
mt     = [   1 1; 1 1];
           % St   Sa  
%data   = [   1.0  0.5  0.5];
data    = [ 1 0.4 0.5 0.0
            1 0.5 0.1 0.5 ];

input   =   struct(   ...
    'numg',            2, ...     % number of groups
    'numm',            1, ...     % number of materials
    'xcm',           xcm, ...     % slab bounds
    'xfm',           xfm, ...     % number of fine meshes
    'ycm',           ycm, ...     % slab bounds
    'yfm',           yfm, ...     % number of fine meshes    
    'mt',             mt, ...     % slab material ids
    'data',         data, ...     % mat comp's
    'src',           src, ...     % volume source
    'ord',             2, ...     % number of ordinates
    'maxit',         100, ...     % max iterations
    'maxerr',       1e-14, ...     % max pointwise phi error
    'adj',             0  ...     % adjoint flag
    );

%[phi,psi,psiV,psiH] = sn_two_d(input);
[KK,psi2,psiH,psiV]=sn_two_d_matrix(input);
[phi,psiC,psiV2,psiH2] = sn_two_d(input);

% check that vertical and horizontal fluxes are the same
max(max(max(max(abs(psiV-psiV2)))))
max(max(max(max(abs(psiH-psiH2)))))

% plot the matrix
spy(KK)
