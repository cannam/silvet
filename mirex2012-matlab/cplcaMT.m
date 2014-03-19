function [w,h,z,u,xa] = cplcaMT( x, K, T, R, w, h, z, u, iter, sw, sh, sz, su, lw, lh, lz, lu, pa)
% function [w,h,xa2] = cplcaMT( x, K, T, R, w, h, z, u, iter, sw, sh, sz, su, lw, lh, lz, lu)
%
% Perform multiple-source, multiple-template SIPLCA for transcription
%
% Inputs:
%  x     input distribution
%  K     number of components
%  T     size of components
%  R     size of sources
%  w     initial value of p(w) [default = random]
%  h     initial value of p(h) [default = random]
%  z     initial value of p(z) [default = random]
%  iter  number of EM iterations [default = 10]
%  sw    sparsity parameter for w [default = 1]
%  sh    sparsity parameter for h [default = 1]
%  sz    sparsity parameter for z [default = 1]
%  lw    flag to update w [default = 1]
%  lh    flag to update h [default = 1]
%  lh    flag to update h [default = 1]
%  pa    source-component activity range [Rx2]
%
% Outputs: 
%  w   p(w) - spectral bases
%  h   p(h) - pitch impulse
%  z   p(z) - mixing matrix for p(h)
%  xa  approximation of input

% Emmanouil Benetos 2011, based on cplca code by Paris Smaragdis


% Sort out the sizes
wc = 2*size(x)-T;
hc = size(x)+T-1;

% Default training iterations
if ~exist( 'iter')
	iter = 10;
end


% Initialize
sumx = sum(x);

if ~exist( 'w') || isempty( w)
    w = cell(R, K);
	for k = 1:K
        for r=1:R
            w{r,k} = rand( T);
            w{r,k} = w{r,k} / sum( w{r,k}(:));
        end
    end
end
if ~exist( 'h') || isempty( h)
    h = cell(1, K);
	for k = 1:K
		h{k} = rand( size(x)-T+1);
		h{k} = h{k};
	end
end
if ~exist( 'z') || isempty( z)
    z = cell(1, K);
	for k = 1:K
		z{k} = rand( size(x,2),1);
		z{k} = z{k};
	end
end
if ~exist( 'u') || isempty( u)
    u = cell(R, K);
	for k = 1:K
        for r=1:R
            if( (pa(r,1) <= k &&  k <= pa(r,2)) )
                u{r,k} = ones( size(x,2),1);
            else
                u{r,k} = zeros( size(x,2),1);
            end
        end;
	end
end

fh = cell(1, K);
fw = cell(R, K);
for k = 1:K
    fh{k} = ones(wc) + 1i*ones(wc);
    for r=1:R
        fw{r,k} = ones(wc) + 1i*ones(wc);
    end;
end;



% Make commands for subsequent multidim operations and initialize fw
fnh = 'c(hc(1)-(T(1)+(1:size(h{k},1))-2),hc(2)-(T(2)+(1:size(h{k},2))-2))';
xai = 'xa(1:size(x,1),1:size(x,2))';
flz = 'xbar(end:-1:1,end:-1:1)';

for k = 1:K
    for r=1:R
        if( (pa(r,1) <= k &&  k <= pa(r,2)) )
            fw{r,k} = fftn( w{r,k}, wc);
        end;
    end;
end;

% Iterate
for it = 1:iter
    
    %disp(['Iteration: ' num2str(it)]);
    
    % E-step
    xa = eps;
    for k = 16:73
        fh{k} = fftn( h{k}, wc);
        for r=1:R
            if( (pa(r,1) <= k &&  k <= pa(r,2)) )
                xa1 = abs( real( ifftn( fw{r,k} .* fh{k})));                
                xa = xa + xa1(1:size(x,1),1:size(x,2)) .*repmat(z{k},1,size(x,1))'.*repmat(u{r,k},1,size(x,1))';
                clear xa1;
            end
        end
    end
    
    xbar = x ./ xa;
    xbar = eval( flz);
    fx = fftn( xbar, wc);
    
    
    % M-step
    for k = 16:73
        
        
        % Update h, z, u
        nh=eps;
        for r=1:R
            if( (pa(r,1) <= k &&  k <= pa(r,2)) )
                c = abs( real( ifftn( fx .* fw{r,k} )));
                nh1 = eval( fnh);
                nh1 = nh1 .*repmat(u{r,k},1,size(h{k},1))';
                nh = nh + nh1;
                
                nhu = eval( fnh);
                nhu = nhu .* h{k};
                nu = sum(nhu)';
                nu = u{r,k} .* nu + eps;
                if lu
                    u{r,k} = nu;
                end;
                
            end;
        end
        nh = h{k} .* (nh.^sh);
        nz = sum(nh)';
        nz = z{k} .* nz + eps;
        
        
        % Assign and normalize
        if lh
            h{k} = nh;
        end
        if lz
            z{k} = nz;
        end
        
        
    end
    
    % Normalize z over t
    if lz
        Z=[]; for i=1:K Z=[Z z{i}]; end;
        Z = Z.^sz;
        Z(1:end,1:15)=0;
        Z(1:end,74:88)=0;
        Z = Z./repmat(sum(Z,2),1,K); z = num2cell(Z,1); %figure; imagesc(imrotate(Z,90));
    end
    
    % Normalize u over z,t
    if lu
        U=[]; for r=1:R U(r,:,:) = cell2mat(u(r,:)); end;
        for i=1:size(U,2) for j=1:size(U,3) U(:,i,j) = U(:,i,j).^su; U(:,i,j) = U(:,i,j) ./ sum(U(:,i,j)); end; end;
        U0 = permute(U,[2 1 3]); u = squeeze(num2cell(U0,1));
    end
    
    % Normalize h over z,t
    H=[]; for k=1:K H(k,:,:) = cell2mat(h(k)); end; H0 = permute(H,[2 1 3]);
    for i=1:size(H0,2) for j=1:size(H0,3) H0(:,i,j) = sumx(j)* (H0(:,i,j) ./ sum(H0(:,i,j))); end; end;
    h = squeeze(num2cell(squeeze(H0),[1 3])); for k=1:K h{k} = squeeze(h{k}); end;
    
    %figure; imagesc(imrotate(xa',90));
    
end

%figure; imagesc(imrotate(xa',90));