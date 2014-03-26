function [w,h,z,xa] = hnmf( x, K, R, iter, sh, sz, w, h, pl, lh, lz)
% function [w,h,z,xa] = hnmf( x, K, R, iter, sh, sz, w, h, pl, lh, lz)
%
% Perform Multi-source NMF
%
% Inputs:
%  x     input distribution
%  K     number of components
%  R     number of sources
%  iter  number of EM iterations [default = 100]
%  sh    sparsity of h
%  sz    sparsity of z
%  w     initial value of w
%  h     initial value of h
%  pl    plot flag
%  lh    update h flag
%  lz    update z flag
%
% Outputs:
%  w   spectral bases
%  h   component activation
%  z   source activation per component
%  xa  approximation of input
%
% Emmanouil Benetos 2011


% Get sizes
[M,N] = size( x);
sumx = sum(x);

% Default training iterations
if ~exist( 'iter')
    iter = 100;
end

% Default plot flag
if ~exist( 'pl')
    pl = 1;
end

% Initialize
if ~exist( 'w') || isempty( w)
    w = rand( M, R, K);
end
for r=1:R
    for k=1:K
        w(:,k,r) = w(:,k,r) ./ sum(w(:,k,r));
    end;
end;
if ~exist( 'h') || isempty( h)
    h = rand( K, N);
end
n=1:N;
h(:,n) = repmat(sumx(n),K,1) .* (h(:,n) ./ repmat( sum( h(:,n), 1), K, 1));
if ~exist( 'z') || isempty( z)
    z = rand( R, K, N);
end
for k=1:K
    for n=1:N
        z(:,k,n) = z(:,k,n) ./ sum(z(:,k,n));
    end;
end;



% Iterate
for it = 1:iter
    
    % E-step
    zh = z .* permute(repmat(h,[1 1 R]),[3 1 2]); %% z is the source activation distribution, h the component (pitch) activation
    xa=eps;
    for r=1:R
        for k=1:K
            xa = xa + w(:,k,r) * squeeze(zh(r,k,:))';
        end; 
    end;
    Q = x ./ xa;
    
    % M-step (update h,z)
    if  (lh && lz)
        nh=zeros(K,N);
        for k=1:K
            for r=1:R
                nh(k,:) = nh(k,:) + squeeze(z(r,k,:))' .* (squeeze(w(:,k,r))' * Q);
                nz = h(k,:) .* (squeeze(w(:,k,r))' * Q);
                nz = nz .* squeeze(z(r,k,:))';
                z(r,k,:) = nz;
            end;
        end;
        nh = h .* nh;
    end
    
    
    % Assign and normalize
    k=1:K;
    n=1:N;
    if lh
        nh = nh.^sh;
        h(:,n) = repmat(sumx(n),K,1) .* (nh(:,n) ./ repmat( sum( nh(:,n), 1), K, 1));
    end
    if lz
        z = z.^sz;
        z(:,k,n) = z(:,k,n) ./ repmat( sum( z(:,k,n), 1), R, 1);
    end
    
end

% Show me
if pl
    subplot(3, 1, 1), imagesc(x), axis xy
    subplot(3, 1, 2), imagesc(xa), axis xy
    subplot(3, 1, 3), imagesc(h), axis xy
end
