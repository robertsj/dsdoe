function [KK,psi,psiH,psiV] = sn_two_d_matrix(in)

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This function solves the 2-D multigroup SN equations given input via a 
% direct linear solve. It's  output is the forward (or adjoint) angular 
% fluxes at the cell edges and the cell-centered scalar flux.  It has been 
% verified against PARTISN for some simple problems.
% ** last modified by J. Roberts, 05/06/2010
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
idxprt = 0;
tic
psi=0;psiH=0;psiV=0;
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
w = 0.25*w;  % weight's sum to unity per octant; normalize so sum(w)=1;
if ord == 44
    ord = 4;
end

N = 0.5*(ord+2)*ord;

Q = zeros(  sum(in.xfm), sum(in.yfm), N, numg );

%bound       =   [ord/2+1 ord; 1 ord/2];
gbound      =   [ 1  numg];
git         =   1;

if in.adj == 1
        git         = -1;
        gbound      = [numg 1];
        %bound       = flipud(bound);
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

nx = sum(in.xfm); % number of x fine meshes
ny = sum(in.yfm); % number of y fine meshes
sigt=zeros(nx,ny,numg);
sigs=zeros(nx,ny,numg,numg);
for g = gbound(1):git:gbound(2)
    for i = 1:nx
        for j = 1:ny
            m = mtt(i,j);
            sigt(i,j,g)   = in.data((m-1)*numg+g,1);
            sigs(i,j,g,:) = in.data((m-1)*numg+g,3:end); 
            % gives sigs( location, groupG' to groupG ) like my excel matrix
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           "K" MATRIX CONSTRUCTION                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% -------------------------------------- DIAGONALS INITIALIZATION

bw = 2*numg*N*(nx+1);                       % matrix bandwidth
kC = ones( (2*nx*ny+nx+ny)*N*numg, 1 );     % central diagonal
kL = zeros( (2*nx*ny+nx+ny)*N*numg, bw-1 ); % lower diagonals
kU = zeros( (2*nx*ny+nx+ny)*N*numg, bw-1 ); % upper diaganals

% make an index vector corresponding diag index to flxgrp and flxord
% the vector is of length bw+, and is in the form
%   | H1,j-1/2,n's H2... H3// |   V(1/2,j,n's) V(3/2).... | H1,j+1/2...|
% ie it holds all possible indices for a given row for a block; below,
% we'll set the first index and go up through bw total indices
flxindex = zeros( N*numg*nx+N*numg*(nx+1), 3);  % 1-->flxord, 2-->flxgrp, 3-->0/1 for H and V (temporary)

% do the H portion
k=0;
for i = 1:N*nx
    for j = 1:numg
        k=k+1;
        ii = mod(i,N); if (ii==0), ii=N; end
        flxindex(k,1)=ii;
        flxindex(k,2)=j;
        flxindex(k,3)=0; % for H
    end
end
for i = 1:N*(nx+1)
    for j = 1:numg
        k=k+1;
        ii = mod(i,N); if (ii==0), ii=N; end
        flxindex(k,1)=ii;
        flxindex(k,2)=j;
        flxindex(k,3)=1; % for V
    end
end

% for the matrix construction, we note the fundamental blocks correspond to 
% rows of the vertical fluxes.  All fluxes for the rows are sandwiched by
% corresponding lower and upper rows of the horizontal fluxes.

% description of indices
%   r = actual row index
%   n = angle, i.e. mu(n) and eta(n)
%   g = group
%   k = column index for upper and lower diagonals
%   kk = maps k to the indices of the fluxes acted upon (via flxindex)
%   nn and gg = angle and group index from flxindex
%   i = the (centered) x coordinate

% -------------------------------------------------------------------------
% -------------------------------------- H BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:ny+1 % loop over all y divisions
for i = 1:nx % loop over all x divisions

r = (i-1)*N*numg + (j-1)*(N*numg*nx+N*numg*(nx+1));

for n = 1:N/4 % for mu < 0, eta < 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for g = 1:numg 
        r = r + 1; % keep track of the actual row index
        if (j < ny + 1)
        % CENTRAL DIAGONAL - non unity values, "minus gamma"
        kC(r) = kC(r) - 0.5*w(n)*sigs(i,j,g,g) / ...
                         (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
        %kC(1)
        % UPPER DIAGONALS
        % This has three components, 1) contiguous H block i.e. those connected
        %   to the central diagonal, 2) the vblock to my right, and 3) the h
        %   block after that v block
        % 1) CONTIGUOUS H BLOCK
        for k = 1:N*numg-((n-1)*numg+g) % length of block past main diagonal        
            kk = numg*(n-1)+g+k;
            nn = flxindex(kk,1); % flxord
            gg = flxindex(kk,2); % flxgrp
            kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                         (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
            if idxprt == 1, kU(r,k) = nn+100*gg; end
        end
        % 2) THE V BLOCK
        % Now k (the column index) must be set to allign with the beginning of 
        %   the verticle block 
        kstart = nx*numg*N-((n-1)*numg+g)+1;
        for k = kstart:kstart+2*N*numg-1 % a full row
            kk = numg*(n-1)+g+k;
            nn = flxindex(kk,1); % flxord
            gg = flxindex(kk,2); % flxgrp
            %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
            kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                         (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
            if idxprt == 1, kU(r,k) = nn+100*gg; end
        end
        % now to subtract the "A" term, following the boxed equations
        kU(r,kstart+(N+n-1)*numg+g-1) = kU(r,kstart+N*numg+g-1) - 4*abs(mu(n))/dx(i) / ...
                          (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
        % 3) THE NEXT H BLOCK
        % Again, reset k to allign with the second H block
        kstart = nx*numg*N + (nx+1)*numg*N -((n-1)*numg+g)+1;
        for k = kstart:kstart + N*numg - 1 % a full row
            kk = k - (nx*numg*N + (nx+1)*numg*N -((n-1)*numg+g));
            nn = flxindex(kk,1); % flxord
            gg = flxindex(kk,2); % flxgrp
            %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
            kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                         (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
            if idxprt == 1, kU(r,k) = nn+100*gg; end
        end  
        % now to add "1-B" term
        kU(r,kstart+(n-1)*numg+g-1) = kU(r,kstart+(n-1)*numg+g-1) + 1 - 4*abs(eta(n))/dy(j) / ...
                          (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
 
        %if (N*numg/4) > 2, LOWER DIAGONALS (only for ord=2,numg=1 are there no)
        if (N*numg/4) > 2 && n*g>1
            % NOTE: lower diagonal 1 is the first below central, 2 is second,
            % etc; in other words, going along a row in reverse!!!
            kmax=(n-1)*numg+g-1;
            for k = 1:kmax % here, we start with flxang 1:<central
                nn = flxindex(k,1);
                gg = flxindex(k,2);
                kL(r,kmax-k+1) = -0.5*w(nn)*sigs(i,j,gg,g) / ...
                         (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                if idxprt == 1, kL(r,kmax-k+1) = nn+100*gg; end
            end
        end
        end
    end
end

for n = N/4+1:N/2 % for mu > 0, eta < 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for g = 1:numg 
        r = r + 1; % keep track of the actual row index
        
        if (j < ny+1)
        % CENTRAL DIAGONAL - non unity values, "minus gamma"
        kC(r) = kC(r) - 0.5*w(n)*sigs(i,j,g,g) / ...
                         (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
        %kC(1)
        % UPPER DIAGONALS
        % This has three components, 1) contiguous H block i.e. those connected
        %   to the central diagonal, 2) the vblock to my right, and 3) the h
        %   block after that v block
        % 1) CONTIGUOUS H BLOCK
        for k = 1:N*numg-((n-1)*numg+g) % length of block past main diagonal        
            kk = numg*(n-1)+g+k;
            nn = flxindex(kk,1); % flxord
            gg = flxindex(kk,2); % flxgrp
            kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                         (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
            if idxprt == 1, kU(r,k) = nn+100*gg; end
        end
        %kU(1:5,1:5)
        % 2) THE V BLOCK
        % Now k (the column index) must be set to allign with the beginning of 
        %   the verticle block
        kstart = nx*numg*N-((n-1)*numg+g)+1;
        for k = kstart:kstart+2*N*numg-1 % a full row
            kk = numg*(n-1)+g+k;
            nn = flxindex(kk,1); % flxord
            gg = flxindex(kk,2); % flxgrp
            flxindex(kk,3);
            %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
            kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                         (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
            if idxprt == 1, kU(r,k) = nn+100*gg; end
        end
        % now to subtract the "A" term, following the boxed equations
        kU(r,kstart+(n-1)*numg+g-1) = kU(r,kstart+(n-1)*numg+g-1) - 4*abs(mu(n))/dx(i) / ...
                          (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
        % 3) THE NEXT H BLOCK
        % Again, reset k to allign with the second H block
        kstart = nx*numg*N + (nx+1)*numg*N -((n-1)*numg+g)+1;
        for k = kstart:kstart + N*numg - 1 % a full row
            kk = k - (nx*numg*N + (nx+1)*numg*N -((n-1)*numg+g));
            nn = flxindex(kk,1); % flxord
            gg = flxindex(kk,2); % flxgrp
            %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
            kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                         (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
            if idxprt == 1, kU(r,k) = nn+100*gg; end
        end  
        % now to add "1-B" term
        kU(r,kstart+(n-1)*numg+g-1) = kU(r,kstart+(n-1)*numg+g-1) + 1 - 4*abs(eta(n))/dy(j) / ...
                          (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
 

        % NOTE: lower diagonal 1 is the first below central, 2 is second,
        % etc; in other words, going along a row in reverse!!!
        kmax=(n-1)*numg+g-1;
        for k = 1:kmax % here, we start with flxang 1:<central
          nn = flxindex(k,1);
          gg = flxindex(k,2);
          kL(r,kmax-k+1) = -0.5*w(nn)*sigs(i,j,gg,g) / ...
              (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
          if idxprt == 1, kL(r,kmax-k+1) = nn+100*gg; end
        end
    
        end
    end
end

for n = N/2+1:3*N/4 % for mu < 0, eta > 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for g = 1:numg
        if j > 1
            r = r + 1; % keep track of the actual row index
            % CENTRAL DIAGONAL - non unity values, "minus gamma"
            
            kC(r) = kC(r) - 0.5*w(n)*sigs(i,j-1,g,g) / ...
                (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
            %kC(1)
            % LOWER
            % This has three components, 1) contiguous H block i.e. those connected
            %   to the central diagonal, 2) the vblock to my right, and 3) the h
            %   block after that v block
            % 1) CONTIGUOUS H BLOCK (n-1)*numg+g-1+numg*N;
            kmax = (n-1)*numg-1+g;%n*numg-1;%+%numg(*N/4;
            for k = 1:kmax % length of block past main diagonal
                nn = flxindex(k,1); % flxord
                gg = flxindex(k,2); % flxgrp
                kL(r,kmax-k+1) =  -0.5*w(nn)*sigs(i,j-1,gg,g) / ...
                             (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
                if idxprt == 1, kL(r,kmax-k+1) = nn+100*gg; end
            end
          
            % 2) THE V BLOCK
            % Now k (the column index) must be set to allign with the beginning of
            %   the verticle block nx*numg*N-((n-1)*numg+g)+1
            kstart = N*numg*(nx)+(n-1)*numg+g - N*numg;
            kmax = kstart+2*N*numg-1;
            for k = kstart:kmax% a full row
                kk = kmax-k+1;%numg*(n-1)+g+k;
                nn = flxindex(kk,1); % flxord
                gg = flxindex(kk,2); % flxgrp
                flxindex(kk,3);
                %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                kL(r,k) =  -0.5*w(nn)*sigs(i,j-1,gg,g) / ...
                    (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
                if idxprt == 1, kL(r,k) = nn+100*gg; end
            end
            % now to subtract the "A" term, following the boxed equations
            kL(r,kstart+(N*numg-((n-1)*numg+g))) = ...
                kL(r,kstart+(N*numg-((n-1)*numg+g))) - 4*abs(mu(n))/dx(i) / ...
                (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
            % 3) THE NEXT H BLOCK
            kstart = N*numg*(nx)+(n-1)*numg+g - N*numg + (nx+1)*N*numg;
            kmax = kstart+N*numg-1;
            for k = kstart:kmax % a full row
                kk = kmax-k+1;% = mod(numg*(n-1)+g+k,N*numg*nx+N*numg*(nx+1));
                %if kk == 0, kk = N*numg*nx+N*numg*(nx+1); end
                nn = flxindex(kk,1); % flxord
                gg = flxindex(kk,2); % flxgrp
                flxindex(kk,3);
                %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                kL(r,k) =  -0.5*w(nn)*sigs(i,j-1,gg,g) / ...
                    (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
                if idxprt == 1, kL(r,k) = nn+100*gg; end
            end
            kL(r,kstart+(N*numg-((n-1)*numg+g))) = ...
                kL(r,kstart+(N*numg-((n-1)*numg+g))) + 1 - 4*abs(eta(n))/dy(j-1) / ...
                (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1)); 
            
            kmax=N*numg-((n-1)*numg+g);
            for k = 1:kmax % here, we start with flxang 1:<central
                kk = numg*(n-1)+g+k;
                nn = flxindex(kk,1);
                gg = flxindex(kk,2);
                kU(r,k) =  -0.5*w(nn)*sigs(i,j-1,gg,g) / ...
                    (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
                if idxprt == 1, kU(r,k) = nn+100*gg; end
            end
        end
    end
end

for n = 3*N/4+1:N % for mu > 0, eta > 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for g = 1:numg
        
            r = r + 1; % keep track of the actual row index
            % CENTRAL DIAGONAL - non unity values, "minus gamma"
         if j > 1   
            kC(r) = kC(r) - 0.5*w(n)*sigs(i,j-1,g,g) / ...
                (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
            %kC(1)
            % LOWER
            % This has three components, 1) contiguous H block i.e. those connected
            %   to the central diagonal, 2) the vblock to my right, and 3) the h
            %   block after that v block
            % 1) CONTIGUOUS H BLOCK (n-1)*numg+g-1+numg*N;
            kmax = (n-1)*numg-1+g;%n*numg-1;%+%numg(*N/4;
            for k = 1:kmax % length of block past main diagonal
                nn = flxindex(k,1); % flxord
                gg = flxindex(k,2); % flxgrp
                kL(r,kmax-k+1) =  -0.5*w(nn)*sigs(i,j-1,gg,g) / ...
                             (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
                if idxprt == 1, kL(r,kmax-k+1) = nn+100*gg; end
            end
            % 2) THE V BLOCK
            % Now k (the column index) must be set to allign with the beginning of
            %   the verticle block
            kstart = N*numg*(nx)+(n-1)*numg+g - N*numg;
            kmax = kstart+2*N*numg-1;
            for k = kstart:kmax% a full row
                kk = kmax-k+1;%numg*(n-1)+g+k;
                nn = flxindex(kk,1); % flxord
                gg = flxindex(kk,2); % flxgrp
                flxindex(kk,3);
                %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                kL(r,k) =  -0.5*w(nn)*sigs(i,j-1,gg,g) / ...
                    (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
                if idxprt == 1, kL(r,k) = nn+100*gg; end
            end
            % now to subtract the "A" term, following the boxed equations
            kL(r,kstart+N*numg+(N*numg-((n-1)*numg+g))) = ...
                     kL(r,kstart+N*numg+(N*numg-((n-1)*numg+g))) - 4*abs(mu(n))/dx(i) / ...
                     (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1)); 
            
            % 3) THE NEXT H BLOCK
            % Again, reset k to allign with the second H block
            kstart = N*numg*(nx)+(n-1)*numg+g - N*numg + (nx+1)*N*numg;
            kmax = kstart+N*numg-1;
            for k = kstart:kmax % a full row
                kk = kmax-k+1;% mod(numg*(n-1)+g+k,N*numg*nx+N*numg*(nx+1));
                %if kk == 0, kk=N*numg*nx+N*numg*(nx+1); end
                nn = flxindex(kk,1); % flxord
                gg = flxindex(kk,2); % flxgrp
                flxindex(kk,3);
                %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                kL(r,k) =  -0.5*w(nn)*sigs(i,j-1,gg,g) / ...
                    (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
                if idxprt == 1, kL(r,k) = nn+100*gg; end
            end
            % now to add "1-B" term
            kL(r,kstart+(N*numg-((n-1)*numg+g))) = ...
                     kL(r,kstart+(N*numg-((n-1)*numg+g))) + 1 - 4*abs(eta(n))/dy(j-1) / ...
                     (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1)); 
            
            %if (N*numg/4) > 2, UPPER DIAGONALS (only for ord=2,numg=1 are
            %there no)
            % NOTE: lower diagonal 1 is the first below central, 2 is second,
            % etc; in other words, going along a row in reverse!!!
            kmax=N*numg-((n-1)*numg+g);
            for k = 1:kmax % here, we start with flxang 1:<central
                kk = numg*(n-1)+g+k;
                nn = flxindex(kk,1);
                gg = flxindex(kk,2);
                kU(r,k) =  -0.5*w(nn)*sigs(i,j-1,gg,g) / ...
                    (sigt(i,j-1,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j-1));
                if idxprt == 1, kU(r,k) = nn+100*gg; end
            end
        end
        
    end
end

end % nx loop
end % ny loop

% -------------------------------------------------------------------------
% -------------------------------------- V BLOCKS %%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:ny
    for i = 1:nx+1
        r = (i-1)*N*numg + (j-1)*(N*numg*nx+N*numg*(nx+1)) + nx*numg*N;
        % same row defininit as above for H but with the 1st h block lengt
        
        for n = 1:N/4 % for mu < 0, eta < 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for g = 1:numg
                r = r + 1; % keep track of the actual row index
                
                if (i < nx + 1)
                    
                    % CENTRAL DIAGONAL - non unity values, "minus gamma"
                    kC(r) = kC(r) - 0.5*w(n)*sigs(i,j,g,g) / ...
                        (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                    
                    % UPPER DIAGONALS
                    % 1) CONTIGUOUS V BLOCK
                    for k = 1:2*N*numg-((n-1)*numg+g)
                        kk = numg*(n-1)+g+k;
                        nn = flxindex(kk,1);
                        gg = flxindex(kk,2);
                        kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                            (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kU(r,k) = nn+100*gg; end
                    end
                    kU(r,numg*N) = 1 + kU(r,numg*N) - 4*abs(mu(n))/dx(i) / ...
                        (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                    
                    % 2) THE UPPER DIAGONAL H BLOCK
                    kstart = (nx+1)*numg*N -((n-1)*numg+g)+1;
                    for k = kstart:kstart+N*numg-1 % a full row
                        kk = k - ( (nx+1)*numg*N -((n-1)*numg+g));
                        nn = flxindex(kk,1); % flxord
                        gg = flxindex(kk,2); % flxgrp
                        %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                        kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                            (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kU(r,k) = nn+100*gg; end
                    end
                    % now to subtract the "B" term, following the boxed equations
                    kU(r,kstart+(n-1)*numg+g-1) = kU(r,kstart+(n-1)*numg+g-1) - 4*abs(eta(n))/dy(j) / ...
                        (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                    
                    % 3) THE LOWER DIAGONAL H BLOCK
                    % Now k (the column index) must be set to allign with the beginning of
                    %   the verticle block
                    kstart = N*numg*(nx)+(n-1)*numg+g - N*numg;
                    kmax = kstart+N*numg-1;
                    for k = kstart:kmax% a full row
                        kk = kmax-k+1;%numg*(n-1)+g+k;
                        nn = flxindex(kk,1); % flxord
                        gg = flxindex(kk,2); % flxgrp
                        flxindex(kk,3);
                        %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                        kL(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                            (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kL(r,k) = nn+100*gg; end
                    end
                    % no B term to subtract
                    
                    % LOWER DIAGONALS
                    if (N*numg/4) > 2 && n*g>1
                        kmax=(n-1)*numg+g-1;
                        for k = 1:kmax % here, we start with flxang 1:<central
                            nn = flxindex(k,1);
                            gg = flxindex(k,2);
                            kL(r,kmax-k+1) = -0.5*w(nn)*sigs(i,j,gg,g) / ...
                                (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                            if idxprt == 1, kL(r,kmax-k+1) = nn+100*gg; end
                        end
                    end
                    
                end % if
                
            end
        end
 
        for n = N/4+1:N/2 % for mu > 0, eta < 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%
            for g = 1:numg
                r = r + 1; % keep track of the actual row index
                if (i > 1)
                    
                    % CENTRAL DIAGONAL - non unity values, "minus gamma"
                    kC(r) = kC(r) - 0.5*w(n)*sigs(i-1,j,g,g) / ...
                        (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                    
                    % LOWER DIAGONALS
                    % 1) CONTIGUOUS BLOCK
                    kmax = (n-1)*numg-1+g+N*numg;%n*numg-1;%+%numg(*N/4;
                    for k = 1:kmax % length of block past main diagonal
                        nn = flxindex(k,1); % flxord
                        gg = flxindex(k,2); % flxgrp
                        kL(r,kmax-k+1) =  -0.5*w(nn)*sigs(i-1,j,gg,g) / ...
                                     (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kL(r,kmax-k+1) = nn+100*gg; end
                    end
                    kL(r,N*numg) = 1 + kL(r,N*numg) - 4*abs(mu(n))/dx(i-1) / ...
                        (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                    
                    
                    % 2) THE UPPER DIAGONAL H BLOCK
                    kstart = (nx+1)*numg*N -((n-1)*numg+g)+1-numg*N;
                    for k = kstart:kstart+N*numg-1 % a full row
                        kk = k - ( (nx+1)*numg*N -((n-1)*numg+g)-numg*N);
                        nn = flxindex(kk,1); % flxord
                        gg = flxindex(kk,2); % flxgrp
                        %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                        kU(r,k) =  -0.5*w(nn)*sigs(i-1,j,gg,g) / ...
                            (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kU(r,k) = nn+100*gg; end
                    end
                    kU(r,kstart+(n-1)*numg+g-1) = ...
                        kU(r,kstart+(n-1)*numg+g-1) - 4*abs(eta(n))/dy(j) / ...
                        (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                    
                    % 3) THE LOWER DIAGONAL H BLOCK
                    % Now k (the column index) must be set to allign with the beginning of
                    %   the verticle block
                    kstart = N*numg*(nx)+(n-1)*numg+g;
                    kmax = kstart+N*numg-1;
                    for k = kstart:kmax% a full row
                        kk = kmax-k+1;%numg*(n-1)+g+k;
                        nn = flxindex(kk,1); % flxord
                        gg = flxindex(kk,2); % flxgrp
                        flxindex(kk,3);
                        %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                        kL(r,k) =  -0.5*w(nn)*sigs(i-1,j,gg,g) / ...
                            (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kL(r,k) = nn+100*gg; end
                    end
                    
                    % Main upper diags
                    kmax=N*numg-((n-1)*numg+g);
                    for k = 1:kmax % here, we start with flxang 1:<central
                        kk = numg*(n-1)+g+k;
                        nn = flxindex(kk,1);
                        gg = flxindex(kk,2);
                        kU(r,k) =  -0.5*w(nn)*sigs(i-1,j,gg,g) / ...
                            (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kU(r,k) = nn+100*gg; end
                    end
                    
                end
            end
        end
        
        for n = N/2+1:3*N/4 % for mu < 0, eta > 0 %%%%%%%%%%%%%%%%%%%%%%%%%
                for g = 1:numg
                    r = r + 1;
                    if (i<nx+1)
                        
                    % CENTRAL DIAGONAL - non unity values, "minus gamma"
                    kC(r) = kC(r) - 0.5*w(n)*sigs(i,j,g,g) / ...
                        (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                    
                    % UPPER DIAGONALS
                    % 1) CONTIGUOUS V BLOCK
                    for k = 1:2*N*numg-((n-1)*numg+g) 
                        kk = numg*(n-1)+g+k;
                        nn = flxindex(kk,1);
                        gg = flxindex(kk,2);
                        kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                            (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kU(r,k) = nn+100*gg; end
                    end
                    kU(r,N*numg) = kU(r,N*numg) + 1 - 4*abs(mu(n))/dx(i) / ...
                        (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                    
                    % 2) THE UPPER DIAGONAL H BLOCK
                    kstart = (nx+1)*numg*N -((n-1)*numg+g)+1;
                    for k = kstart:kstart+N*numg-1 % a full row
                        kk = k - ( (nx+1)*numg*N -((n-1)*numg+g));
                        nn = flxindex(kk,1); % flxord
                        gg = flxindex(kk,2); % flxgrp
                        %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                        kU(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                            (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kU(r,k) = nn+100*gg; end
                    end
                    % no "B" term

                    % 3) THE LOWER DIAGONAL H BLOCK
                    % Now k (the column index) must be set to allign with the beginning of
                    %   the verticle block
                    kstart = N*numg*(nx)+(n-1)*numg+g - N*numg;
                    kmax = kstart+N*numg-1;
                    for k = kstart:kmax% a full row
                        kk = kmax-k+1;%numg*(n-1)+g+k;
                        nn = flxindex(kk,1); % flxord
                        gg = flxindex(kk,2); % flxgrp
                        flxindex(kk,3);
                        %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                        kL(r,k) =  -0.5*w(nn)*sigs(i,j,gg,g) / ...
                            (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kL(r,k) = nn+100*gg; end
                    end
                    kL(r,kstart+(N*numg-((n-1)*numg+g))) = ...
                        kL(r,kstart+(N*numg-((n-1)*numg+g)))  - 4*abs(eta(n))/dy(j) / ...
                        (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j)); 
            
                    % LOWER DIAGONALS
                    kmax=(n-1)*numg+g-1;
                    for k = 1:kmax % here, we start with flxang 1:<central
                        nn = flxindex(k,1);
                        gg = flxindex(k,2);
                        kL(r,kmax-k+1) = -0.5*w(nn)*sigs(i,j,gg,g) / ...
                            (sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kL(r,kmax-k+1) = nn+100*gg; end
                    end
                    
                    end    
                end
        end

        for n = 3*N/4+1:N % for mu > 0, eta > 0 %%%%%%%%%%%%%%%%%%%%%%%%%%%
            for g = 1:numg
                r = r + 1;
                if (i > 1)
                    % CENTRAL DIAGONAL - non unity values, "minus gamma"
                    kC(r) = kC(r) - 0.5*w(n)*sigs(i-1,j,g,g) / ...
                        (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                    
                    % LOWER DIAGONALSk
                    % 1) CONTIGUOUS BLOCK
                    kmax = (n-1)*numg-1+g+N*numg;%n*numg-1;%+%numg(*N/4;
                    for k = 1:kmax % length of block past main diagonal
                        nn = flxindex(k,1); % flxord
                        gg = flxindex(k,2); % flxgrp
                        kL(r,kmax-k+1) =  -0.5*w(nn)*sigs(i-1,j,gg,g) / ...
                                     (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kL(r,kmax-k+1) = nn+100*gg; end
                    end
                    kL(r,N*numg) = 1 + kL(r,N*numg) - 4*abs(mu(n))/dx(i-1) / ...
                        (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                     
                    % 2) THE UPPER DIAGONAL H BLOCK
                    kstart = (nx+1)*numg*N -((n-1)*numg+g)+1-numg*N;
                    for k = kstart:kstart+N*numg-1 % a full row
                        kk = k - ( (nx+1)*numg*N -((n-1)*numg+g)-numg*N);
                        nn = flxindex(kk,1); % flxord
                        gg = flxindex(kk,2); % flxgrp
                        %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                        kU(r,k) =  -0.5*w(nn)*sigs(i-1,j,gg,g) / ...
                            (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kU(r,k) = nn+100*gg; end
                    end
                    
                    % 3) THE LOWER DIAGONAL H BLOCK
                    % Now k (the column index) must be set to allign with the beginning of
                    %   the verticle block
                    kstart = N*numg*(nx)+(n-1)*numg+g;
                    kmax = kstart+N*numg-1;
                    for k = kstart:kmax% a full row
                        %kk = numg*(n-1)+g+k;
                        kk = kmax-k+1;
                        nn = flxindex(kk,1); % flxord
                        gg = flxindex(kk,2); % flxgrp
                        flxindex(kk,3);
                        %disp(['kk = ',num2str(kk),' nn = ',num2str(nn),' gg = ',num2str(gg)])
                        kL(r,k) =  -0.5*w(nn)*sigs(i-1,j,gg,g) / ...
                            (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kL(r,k) = nn+100*gg; end
                    end
                    kL(r,kstart+(N*numg-((n-1)*numg+g))) = ...
                     kL(r,kstart+(N*numg-((n-1)*numg+g))) - 4*abs(eta(n))/dy(j) / ...
                     (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j)); 
            
                    % Main upper diags
                    kmax=N*numg-((n-1)*numg+g);
                    for k = 1:kmax % here, we start with flxang 1:<central
                        kk = numg*(n-1)+g+k;
                        nn = flxindex(kk,1);
                        gg = flxindex(kk,2);
                        kU(r,k) =  -0.5*w(nn)*sigs(i-1,j,gg,g) / ...
                            (sigt(i-1,j,g)+2*abs(mu(n))/dx(i-1)+2*abs(eta(n))/dy(j));
                        if idxprt == 1, kU(r,k) = nn+100*gg; end
                    end
                end
            end
        end
        
    end
end


% -------------------------------------------------------------------------
% -------------------------------------- RE-ORDER THE DIAGONALS
L = length(kU(:,1));
% pad the upper diagonals for proper sparse K construction
for i=1:bw-1
   kU( (i+1):end, i) = kU( 1:(L-i), i);
   kU( 1:i, i) = 0.0;
end
% same for lower
for i=1:bw-1
   kL( 1:(L-i), i) = kL( (i+1):end, i);
   %kL( (L-i+1):end, i) = 0.0;
end

% -------------------------------------------------------------------------
% -------------------------------------- BUILD THE SPARSE MATRIX
KK = spdiags([ kL(:,end:-1:1) kC kU ], ...
                    [-(bw-1):(bw-1)], length(kC), length(kC));

% -------------------------------------------------------------------------
% -------------------------------------- BUILD THE RIGHT HAND SIDE

RHS = zeros( (2*nx*ny+nx+ny)*N*numg, 1 );

r = 0;
% build bottom
for i = 1:nx
    for n = 1:N/2
        for g = 1:numg
            r = r + 1;
            RHS(r) = Q(i,1,n,g)*...
                2/(sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
        end
    end
    for n = N/2+1:N
        for g = 1:numg
            r = r + 1;
            % bc's if desired
        end
    end
end

% build middle 
for j = 1:ny
    % vert left hand side
    for n = 1:N/4
        for g = 1:numg
            r = r + 1;
            RHS(r) = Q(1,j,n,g)*...
                2/(sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
        end
    end
    for n = N/4+1:N/2
        for g = 1:numg
            r = r + 1;
            % bc's if desired
        end
    end
    for n = N/2+1:3*N/4
        for g = 1:numg
            r = r + 1;
            RHS(r) = Q(1,j,n,g)*...
                2/(sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
        end
    end
    for n = 3*N/4+1:N
        for g = 1:numg
            r = r + 1;
            % bc's if desired
        end
    end
    % vert middle
    for i = 2:nx
        for n = 1:N
            for g = 1:numg
                r = r + 1;
                RHS(r) = Q(i,j,n,g)*...
                    2/(sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
            end
        end
    end

    % vert right hand side
    for n = 1:N/4
        for g = 1:numg
            r = r + 1;
            % bc times
        end
    end
    for n = N/4+1:N/2 
        for g = 1:numg
            r = r + 1;
            RHS(r) = Q(nx,j,n,g)*...
                2/(sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
        end
    end
    for n = N/2+1:3*N/4
        for g = 1:numg
            r = r + 1;
            % bc
        end
    end
    for n = 3*N/4+1:N
        for g = 1:numg
            r = r + 1;
            RHS(r) = Q(nx,j,n,g)*...
                2/(sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
        end
    end
    if j < ny
       for i = 1:nx
          for n = 1:N
              for g = 1:numg
                r = r + 1;
                RHS(r) = Q(i,j,n,g)*...
                2/(sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j));
              end
          end
       end
    end
end
% build top
for i = 1:nx
    for n = 1:N/2
        for g = 1:numg
            r = r + 1;
            %  bc
        end
    end
    for n = N/2+1:N
        for g = 1:numg
            r = r + 1;
            RHS(r) = Q(i,ny,n,g)*...
                2/(sigt(i,j,g)+2*abs(mu(n))/dx(i)+2*abs(eta(n))/dy(j)); 
        end
    end
end

%for i = 1:length(RHS)
%    disp([num2str(i),'  ',num2str(RHS(i))])
%end
if idxprt == 1
    lala=1;
end
psi = KK\RHS;

%  poo = full(KK);poo=poo(2:2:end,2:2:end); 
%  KK=poo; psi=RHS(2:2:end);
% return

%psi = psi(1:2:end); numg=1;
% get psiH
psiH = zeros(nx,ny+1,N,numg);
for j = 1:ny+1
    r = (j-1)*(2*nx+1)*numg*N;
    for i = 1:nx
        for n = 1:N
            for g = 1:numg
                r = r + 1;
                psiH(i,j,n,g) = psi(r);
            end
        end
    end
end
% get psiV
psiV = zeros(nx+1,ny,N,numg);
for j = 1:ny
    r = nx*numg*N+(j-1)*(2*nx+1)*numg*N;
    for i = 1:nx+1
        for n = 1:N
            for g = 1:numg
                r = r + 1;
                psiV(i,j,n,g) = psi(r);
            end
        end
    end
end


end





