function [phi,psiC,psiV,psiH] = sn_two_d(in)

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This function solves the 2-D multigroup SN equations given input.  It's !
% output is the forward (or adjoint) angular fluxes at the cell edges and !
% the cell-centered scalar flux.  It has been verified against PARTISN for!
% some simple problems.
% ** last modified by J. Roberts, 05/03/2010
%    -- made a function for weights and reordered sweep direction order
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
tic
% input contains
%     numg   = number of energy groups
%     numm   = number of materials
%     xcm    = x coarse divisions
%     xfm    = x fine mesh interval for coarse divisions
%     ycm    = y coarse divisions
%     yfm    = y fine mesh interval for coarse divisions
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
% working variables
%     dx     = fine mesh divisions
%     mtt    = material assignment for fine mesh cells
%     n      = number of fine mesh cells

numg = in.numg;
ord  = in.ord;
[mu,eta,w]=S_2D(ord);
if ord == 44
    ord = 4;
end
N = ord;

numord = 0.5*(N+2)*N;
negflx=0;

Q = zeros(  sum(in.xfm), sum(in.yfm), numord, numg );

bound       =   [ord/2+1 ord; 1 ord/2];
gbound      =   [ 1  numg];
git         =   1;

if in.adj == 1
        git         = -1;
        gbound      = [numg 1];
        bound       = flipud(bound);
end

% -------------------------------------------------------------------------
% -------------------------------------- Discretizations

k = 0;
kk = 0;

dx = zeros(sum(in.xfm),1);
dy = zeros(sum(in.yfm),1);
mtt = zeros(sum(in.xfm),sum(in.yfm));


for i = 1:length(in.xfm)
    dx( (k+1):(k+in.xfm(i))   )  = (in.xcm(i+1)-in.xcm(i))/in.xfm(i);
    for j = 1:length(in.yfm)
        dy( (kk+1):(kk+in.yfm(j))   )  = ...
            (in.ycm(j+1)-in.ycm(j))/in.yfm(j);
        for g=gbound(1):git:gbound(2)
            Q( (k+1):(k+in.xfm(i)), (kk+1):(kk+in.yfm(j)), :, g)  = ...
                in.src(g,i,j);
        end
        mtt( (k+1):(k+in.xfm(i)), (kk+1):(kk+in.yfm(j))   )  = ...
            in.mt(i,j);  % assign mat to each f mesh
        kk = sum(in.yfm(1:j));
    end
    kk = 0;
    k = sum(in.xfm(1:i));
end

nx = sum(in.xfm);
ny = sum(in.yfm);

% -------------------------------------------------------------------------
% -------------------------------------- Transport Coefficients
con1 = zeros(nx,ny,numord,numg);
con2 = zeros(nx,ny,numord,numg);
con3 = zeros(nx,ny,numord,numg);
for g = gbound(1):git:gbound(2)
    for i = 1:nx
        for j = 1:ny
            m = mtt(i,j);
            % constant 1 :
            sig = in.data((m-1)*numg+g,1); % total cross-section
            for n = 1:numord
                con2(i,j,n,g) = 2*abs(mu(n))/dx(i);
                con3(i,j,n,g) = 2*abs(eta(n))/dy(j);
                con1(i,j,n,g) = 1 / (sig + con2(i,j,n,g) + con3(i,j,n,g));
            end
        end
    end
end
      
% initialize psi matrices (could cut down on storage)

psiH = zeros(nx,ny+1,numord,numg);  % horizontal edge flux
psiV = zeros(nx+1,ny,numord,numg); % vertical edge flux
psiC = zeros(nx,ny,numord,numg);
phi  = zeros(nx,ny,numg);

s = Q;

% ----------------- Solution Algorithm ------------------------------------

for g = gbound(1):git:gbound(2)
    
    % ----------------- Convergence Parameters ----------------------------
    phierr = 1; % reset phierr
    it = 0;
    A = 2.0;
    while it < in.maxit && phierr > in.maxerr
        % mu > 0, eta > 0
        for n = 3*numord/4+1:numord
            for i = 2:nx+1           % left to right
                for j = 2:ny+1       % bottom to top
                    psiC(i-1,j-1,n,g) = con1(i-1,j-1,n,g) * ...
                                      (con2(i-1,j-1,n,g)*psiV(i-1,j-1,n,g) + ...
                                       con3(i-1,j-1,n,g)*psiH(i-1,j-1,n,g) + ...
                                       s(i-1,j-1,n,g));
                    if negflx == 1               
                        if psiC(i-1,j-1,n,g)<0, psiC(i-1,j-1,n,g)=0; end    
                    end
                    psiV(i,j-1,n,g) = A*psiC(i-1,j-1,n,g) - psiV(i-1,j-1,n,g);
                    psiH(i-1,j,n,g) = A*psiC(i-1,j-1,n,g) - psiH(i-1,j-1,n,g);
                    if negflx == 1
                    if psiV(i,j-1,n,g)<0, psiV(i,j-1,n,g)=0; end 
                    if psiH(i-1,j,n,g)<0, psiH(i-1,j,n,g)=0; end 
                    end
                   % psiC(i-1,j-1,n,g),psiV(i,j-1,n,g),psiH(i-1,j,n,g)
                end
            end
        end
        % mu < 0, eta > 0
        for n = numord/2+1:3*numord/4
            for i = nx:-1:1      % right to left
                for j = 2:ny+1   % bottom to top
                    psiC(i,j-1,n,g) = con1(i,j-1,n,g) * ...
                                    (con2(i,j-1,n,g)*psiV(i+1,j-1,n,g) + ...
                                     con3(i,j-1,n,g)*psiH(i,j-1,n,g) + ...
                                     s(i,j-1,n,g));
                    if negflx == 1
                    if psiC(i,j-1,n,g)<0, psiC(i,j-1,n,g)=0; end 
                    end            
                    psiV(i,j-1,n,g) = A*psiC(i,j-1,n,g) - psiV(i+1,j-1,n,g);
                    psiH(i,j,n,g)   = A*psiC(i,j-1,n,g) - psiH(i,j-1,n,g);
                    if negflx == 1
                    if psiV(i,j-1,n,g)<0, psiV(i,j-1,n,g)=0; end 
                    if psiH(i,j,n,g)<0, psiH(i,j,n,g)=0; end 
                    end
                    
                end
            end
        end
        % mu > 0, eta < 0
        for n = numord/4+1:numord/2
            for i = 2:nx+1       % left to right
                for j = ny:-1:1  % top to bottom
                    psiC(i-1,j,n,g) = con1(i-1,j,n,g) * ...
                                    (con2(i-1,j,n,g)*psiV(i-1,j,n,g) + ...
                                     con3(i-1,j,n,g)*psiH(i-1,j+1,n,g) + ...
                                     s(i-1,j,n,g)); 
                    if negflx == 1             
                    if psiC(i-1,j,n,g)<0, psiC(i-1,j,n,g)=0; end 
                    end                  
                    psiV(i,j,n,g) = A*psiC(i-1,j,n,g) - psiV(i-1,j,n,g);
                    psiH(i-1,j,n,g) = A*psiC(i-1,j,n,g) - psiH(i-1,j+1,n,g);
                    if negflx == 1
                    if psiV(i,j,n,g)<0, psiV(i,j,n,g)=0; end 
                    if  psiH(i-1,j,n,g)<0, psiH(i-1,j,n,g)=0; end  
                    end             
                end
            end
        end
        % mu < 0, eta < 0
        for n = 1:numord/4
            for i = nx:-1:1      % right to left
                for j = ny:-1:1  % top to bottom
                    psiC(i,j,n,g) = con1(i,j,n,g) * ...
                                  (con2(i,j,n,g)*psiV(i+1,j,n,g) + ...
                                   con3(i,j,n,g)*psiH(i,j+1,n,g) + ...
                                   s(i,j,n,g));        
                    if negflx == 1           
                    if psiC(i,j,n,g)<0, psiC(i,j,n,g)=0; end        
                    end
                    psiV(i,j,n,g) = A*psiC(i,j,n,g) - psiV(i+1,j,n,g);
                    psiH(i,j,n,g) = A*psiC(i,j,n,g) - psiH(i,j+1,n,g);
                    if negflx == 1
                    if psiV(i,j,n,g)<0, psiV(i,j,n,g)=0; end 
                    if psiH(i,j,n,g)<0, psiH(i,j,n,g)=0; end 
                    end
                end
            end
        end    
        % ---- Scalar Flux (L&M eq 4-17, l=0)
        phi0(:,:) = phi(:,:,g);
        phi(:,:,g) = zeros(nx,ny,1);
        for n = 1:numord
            phi(:,:,g) = phi(:,:,g) + 0.25*psiC(:,:,n,g)*w(n);
        end
        % ---- Updated Source Term
        for z   =   g:git:gbound(2) % only down scattering
            for i = 1:nx
                for j = 1:ny
                    s(i,j,:,z) = Q(i,j,:,z); % reset
                end
            end
            if in.adj == 0 % FORWARD TRANSPORT
                for i = 1:nx
                    for j = 1:ny
                        m = mtt(i,j);
                        for gg =  gbound(1):git:gbound(2) % group gg to group z
                            s(i,j,:,z) = s(i,j,:,z) + ...
                                in.data((m-1)*numg+gg,2+z)*phi(i,j,gg);
                        end
                    end
                end
            else % ADJOINT TRANSPORT
                for i = 1:nx
                    for j = 1:ny
                        m = mtt(i,j);
                        for gg = gbound(1):git:gbound(2) % group gg to group z
                            s(i,j,:,z) = s(i,j,:,z) + ...
                                in.data((m-1)*numg+z,2+gg)*phi(i,j,gg);
                        end
                    end
                end
            end
        end
        % ---- Scalar Flux Error Between Iterations
        phierr = max(max(abs(phi(:,:,g)-phi0)./phi(:,:,g)));
        it = it+1;
    end % while
    disp([' Group ', num2str(g),' Iterations:   ',num2str(it)])
end

disp([' Elapsed time: ',num2str(toc)])
disp([' Iterations:   ',num2str(it)])


