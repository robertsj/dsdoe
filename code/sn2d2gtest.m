% This input file creates Figure 4 in the report
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
clear, clc
format short

xcm    = [ 0   2   4    6  8    10 ];
xfm    = [  10   10  10  10   10   ];
ycm    = [ 0   2   4    6   8   10];
yfm    = [  10   10  10  10   10   ];
src(1,:,:)    = [   1   0   0   0    0  
                    0   0   0   0    0  
                    0   0   0   0    0 
                    0   0   0   0    0 
                    0   0   0   0    0 ] ;   
src(2,:,:)    = [   0   0   0   0    0  
                    0   0   0   0    0  
                    0   0   0   0    0 
                    0   0   0   0    0 
                    0   0   0   0    1 ] ;                   
mt     = [   1   1   1   1   1
             1   1   1   1   1 
             1   1   1   1   1 
             1   1   1   1   1 
             1   1   1   1   1 ];
           % St   Sa  
data   = [   1.0  0.5  0.3  0.2000
             1.0  0.6  0.0  0.400];

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
    'ord',            44, ...     % number of ordinates
    'maxit',         100, ...     % max iterations
    'maxerr',       1e-9, ...     % max pointwise phi error
    'adj',             0  ...     % adjoint flag
    );

[phi,psi] = sn_two_d(input);

x = linspace(1,10,50);
[X,Y]=meshgrid(x,x);

subplot(1,2,1)
contourf(X,Y,phi(:,:,1),100)
shading flat
colorbar
%square axis
xlabel('x [cm]'),ylabel('y [cm]'),title('\phi_1(x,y) [n/cm^2-s]')
subplot(1,2,2)
contourf(X,Y,phi(:,:,2),100)
xlabel('x [cm]'),ylabel('y [cm]'),title('\phi_2(x,y) [n/cm^2-s]')
shading flat
colorbar
%square axis