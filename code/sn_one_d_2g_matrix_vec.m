function [KK,Q] = sn_one_d_2g_matrix_vec(in)
%profile on
%tic
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This function casts the 1-D multigroup SN equations into matrix form    !
% for direct solution via elimination or krylov solvers                   !
% ** last modified by J. Roberts, 5/13/2010                               !
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% TO DO:
%   Modify numerical constants to minimize divisions 
%   add computation of the scalar flux
%   add reflective conditions
%   add incident conditions

numg = in.numg;
N    = in.ord; 

gbound      =   [ 1  numg];
git         =   1;

% -------------------------------------------------------------------------
% -------------------------------------- Discretizations
[mu,w] = S_1D(N); % get cosines and weights
w = 0.5*w;        % normalize weights to unity

nfm = sum(in.xfm); % number of fine meshes

ssrc = zeros(nfm,N,numg);

j = 0;
for i = 1:length(in.xfm)
    dx( (j+1):(j+in.xfm(i))   )  = (in.xcm(i+1)-in.xcm(i))/in.xfm(i);
    for g=gbound(1):git:gbound(2)
        ssrc( (j+1):(j+in.xfm(i)), :, g)  = in.src(g,i);
    end
    mtt( (j+1):(j+in.xfm(i))   )  = in.mt(i);  % assign mat to each f mesh
    j = sum(in.xfm(1:i));
end

sigt = zeros(nfm,numg);
sigs = zeros(nfm,numg,numg);
for g = gbound(1):git:gbound(2)
    for k = 1:nfm
        m = mtt(k);
        sigt(k,g)   = in.data((m-1)*numg+g,1);
        sigs(k,g,:) = in.data((m-1)*numg+g,3:end); 
    end
end

%time(1) = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           "K" MATRIX CONSTRUCTION                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% -------------------------------------- DIAGONALS INITIALIZATION
bw = 2*numg*N; % matrix bandwidth
kC = ones( (nfm+1)*N*numg, 1 );     % central diagonal
kL = zeros( (nfm+1)*N*numg, bw-1 ); % lower diagonals
kU = zeros( (nfm+1)*N*numg, bw-1 ); % upper diaganals

% make an index vector corresponding diag index to flxgrp and flxord
flxindex=zeros(2*N*numg,2);
k=0;
for i = 1:2*N
    for j = 1:numg
        k=k+1;
        ii = mod(i,N); if (ii==0), ii=N; end
        flxindex(k,1)=ii;
        flxindex(k,2)=j;
    end
end

%time(2) = toc;

% -------------------------------------------------------------------------
% -------------------------------------- FORWARD BLOCK, FIRST REGION
j = 0;
for n = 1:N/2 % for mu < 0
    for g = 1:numg 
        j = j + 1;
        % CENTRAL DIAGONAL
        kC(j) = kC(j) - 0.5*sigs(1,g,g)*w(n)*dx(1)/...
                           abs(mu(n))/(1.0+0.5*dx(1)*sigt(1,g)/abs(mu(n)));

        for k = 1:2*N*numg-((n-1)*numg+g) 
            
            kk = numg*(n-1)+g+k;
            ii = flxindex(kk,1);
            gg = flxindex(kk,2);

            kU(j,k) = -0.5*sigs(1,gg,g)*w(ii)*dx(1)/...
                abs(mu(n))/(1.0+0.5*dx(1)*sigt(1,g)/abs(mu(n)));

        end

        % update with the "alpha" term
        kU(j,N*numg) = kU(j,N*numg)-(1-sigt(1,g)*dx(1)/2/abs(mu(n)))/...
            (1+sigt(1,g)*dx(1)/2/abs(mu(n)));
        
        % if N*numg > 2, LOWER DIAGONALS
        if (N*numg > 2 && n*g > 1)
            kmax=(n-1)*numg+g-1;
            for k = 1:kmax % here, we start with flxang 1:<central
                ii = flxindex(k,1);
                gg = flxindex(k,2);
                kL(j,kmax-k+1) = -0.5*sigs(1,gg,g)*w(ii)*dx(1)/abs(mu(n))/...
                    (1.0+0.5*dx(1)*sigt(1,g)/abs(mu(n)));
            end
        end
        
    end
end
%time(3) = toc;

% -------------------------------------------------------------------------
% -------------------------------------- BACKWARD BLOCK, LAST REGION
j = nfm*(N*numg)+N/2*numg; % starts it in the last block
for n = N/2+1:N
    for g = 1:numg
        
        j = j+1; 
        
        % CENTRAL DIAGONAL - non unity values
        kC(j) = kC(j) - 0.5*sigs(end,g,g)*w(n)*dx(end)/abs(mu(n))...
            /(1.0+0.5*dx(end)*sigt(end,g)/abs(mu(n)));
        
        % LOWER DIAGONAL
        kmax=(n-1)*numg+g-1+numg*N;
        for k = 1:kmax
            ii = flxindex(k,1);
            gg = flxindex(k,2);
            kL(j,kmax-k+1) = -0.5*sigs(end,gg,g)*w(ii)*dx(end)/abs(mu(n))/...
                (1.0+0.5*dx(end)*sigt(end,g)/abs(mu(n)));
        end
        kL(j,N*numg) = kL(j,N*numg) - (1-sigt(end,g)*dx(end)/2/abs(mu(n)))/...
            (1+sigt(end,g)*dx(end)/2/abs(mu(n)));
        
        % if N > 2, UPPER DIAGONALS
        if (N*numg > 2 && n*g > 1)
            for k = 1:2*N*numg-((n-1)*numg+g) 
                
                kk = numg*(n-1)+g+k;
                ii = flxindex(kk,1);
                gg = flxindex(kk,2);

                kU(j,k) = -0.5*sigs(end,gg,g)*w(ii)*dx(end)/...
                    abs(mu(n))/(1.0+0.5*dx(end)*sigt(end,g)/abs(mu(n)));

            end
        end
  
    end
end
%time(4) = toc;
% -------------------------------------------------------------------------
% -------------------------------------- REST OF THE BLOCKS
for m = 2:nfm % all but the first and last regions
    
    % THE FORWARD BLOCK
    j = (m-1)*N*numg;
    for n = 1:N/2
        for g = 1:numg
            j=j+1;
            
            % CENTRAL DIAGONAL - non unity values
            kC(j) = kC(j) - 0.5*sigs(m,g,g)*w(n)*dx(m)/abs(mu(n))/...
                (1.0+0.5*dx(m)*sigt(m,g)/abs(mu(n)));
            
            % UPPER DIAGONALS
            for k = 1:2*N*numg-((n-1)*numg+g)
                kk = numg*(n-1)+g+k;
                ii = flxindex(kk,1);
                gg = flxindex(kk,2);
                kU(j,k) = -0.5*sigs(1,gg,g)*w(ii)*dx(1)/...
                    abs(mu(n))/(1.0+0.5*dx(m)*sigt(m,g)/abs(mu(n)));
            end
            % update with the "alpha" term
            kU(j,N*numg) = kU(j,N*numg)-(1-sigt(m,g)*dx(m)/2/abs(mu(n)))/...
                (1+sigt(m,g)*dx(m)/2/abs(mu(n)));
            
            % if N*numg > 2, LOWER DIAGONALS
            if (N*numg > 2 && n*g > 1)
                kmax=(n-1)*numg+g-1;
                for k = 1:kmax % here, we start with flxang 1:<central
                    ii = flxindex(k,1);
                    gg = flxindex(k,2);
                    kL(j,kmax-k+1) = -0.5*sigs(m,gg,g)*w(ii)*dx(m)/...
                        abs(mu(n))/(1.0+0.5*dx(m)*sigt(m,g)/abs(mu(n)));
                end
            end
            
        end
    end
    
    % THE BACKWARD BLOCK
    j = (m-1)*N*numg + N/2*numg;
    for n = N/2+1:N
        for g = 1:numg
            j=j+1;
            
            % CENTRAL DIAGONAL - non unity values
            kC(j) = kC(j) - 0.5*sigs(m-1,g,g)*w(n)*dx(m-1)/abs(mu(n))...
                /(1.0+0.5*dx(m-1)*sigt(m-1,g)/abs(mu(n)));
            
            % LOWER DIAGONAL
            kmax=(n-1)*numg+g-1+numg*N;
            for k = 1:kmax
                ii = flxindex(k,1);
                gg = flxindex(k,2);
                kL(j,kmax-k+1) = -0.5*sigs(m-1,gg,g)*w(ii)*dx(m-1)/abs(mu(n))/...
                    (1.0+0.5*dx(m-1)*sigt(m-1,g)/abs(mu(n)));
            end
            kL(j,N*numg) = kL(j,N*numg) - (1-sigt(m-1,g)*dx(m-1)/2/abs(mu(n)))/...
                (1+sigt(m-1,g)*dx(m-1)/2/abs(mu(n)));
            
            %if N > 2, UPPER DIAGONALS
            if (N*numg > 2 && n*g > 1)
                for k = 1:2*N*numg-((n-1)*numg+g+N*numg)
                    kk = numg*(n-1)+g+k;
                    ii = flxindex(kk,1);
                    gg = flxindex(kk,2);
                    kU(j,k) = -0.5*sigs(m-1,gg,g)*w(ii)*dx(m-1)/...
                        abs(mu(n))/(1.0+0.5*dx(m-1)*sigt(m-1,g)/abs(mu(n)));
                end
            end
        
        end
    end    
end % regions
%time(5) = toc;
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
   kL( (L-i+1):end, i) = 0.0;
end
%time(6) = toc;
% -------------------------------------------------------------------------
% -------------------------------------- BUILD THE SPARSE MATRIX
KK = spdiags([ kL(:,end:-1:1) kC kU ], ...
                    -(bw-1):(bw-1), length(kC), length(kC));
%time(7) = toc;
% -------------------------------------------------------------------------
% -------------------------------------- SOURCE TERM (i.e. RHS)
Q = zeros( numg*N*(nfm+1), 1 ); % Q = beta(mat,n)*S(mat,n);
% first edge
j=0;
for i = 1:N/2
    for g = 1:numg
        j=j+1;
        Q(j) = dx(1)/abs(mu(i))/(1+0.5*dx(1)*sigt(1,g)/abs(mu(i))) * ssrc(1,g);
    end
end
% last edge
j = nfm*N*numg+N/2*numg;
for i = N/2+1:N
    for g = 1:numg
        j=j+1;
        Q(j) = dx(end)/abs(mu(i))/(1+0.5*dx(end)*sigt(end,g)/abs(mu(i))) * ssrc(end,g);
    end
end
% other edges
for m = 2:nfm
    j=(m-1)*N*numg;
    for i = 1:N/2
        for g = 1:numg
            j=j+1;
            Q(j) =  dx(m)/abs(mu(i))/(1+0.5*dx(m)*sigt(m,g)/abs(mu(i))) * ssrc(m,g);
        end
    end
    j=(m-1)*N*numg+N/2*numg;
    for i = N/2+1:N
        for g = 1:numg
            j=j+1;
            Q(j) =  dx(m-1)/abs(mu(i))/(1+0.5*dx(m-1)*sigt(m-1,g)/abs(mu(i))) * ssrc(m-1,g);
        end
    end
    
end
%time(8) = toc;

% for i = 1:8
%     if i==1
%         disp(['%time ',num2str(i),' = ',num2str(t(i))])
%     else
%         disp(['%time ',num2str(i),' = ',num2str(t(i)-t(i-1))])
%     end
%     
% end

%profile report
end