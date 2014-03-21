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


%% for the transcription application,
%% x -> noise-reduced constant Q. In the application this is a 2-sec,
%%   100-col segment with 2 zeros at top and bottom, so 549x100
%% K -> 88, number of notes
%% T -> [545 1], a two-element array: 545 is the length of each
%%   template, but why 1?
%% R -> 10, number of instruments
%% w -> a 10x88 cell array, in which w{instrument,note} is a 545x1
%%   array containing the template for the given instrument and note
%%   number
%% h -> empty
%% z -> empty
%% u -> empty
%% iter -> a parameter for the program, 12 in the mirex submission
%% sw -> 1
%% sh -> 1
%% sz -> 1.1
%% su -> 1.3, not documented above, presumably sparsity for u
%% lw -> 0, don't update w
%% lh -> 1, do update h
%% lz -> 1, do update z
%% lu -> 1, not documented above, presumably do update u
%% pa -> a 10x2 array in which pa(instrument,1) is the lowest note
%%   expected for that instrument and pa(instrument,2) is the highest


% Sort out the sizes

wc = 2*size(x)-T; %% works out to 553x199
hc = size(x)+T-1; %% works out to 1093x100

% Default training iterations
if ~exist( 'iter')
	iter = 10;
end


% Initialize
sumx = sum(x); %% for later normalisation

if ~exist( 'w') || isempty( w)
    %% doesn't happen, w was provided (it's the template data)
    w = cell(R, K);
	for k = 1:K
        for r=1:R
            w{r,k} = rand( T);
            w{r,k} = w{r,k} / sum( w{r,k}(:));
        end
    end
end
if ~exist( 'h') || isempty( h)
    %% does happen, h was not provided
    h = cell(1, K);
	for k = 1:K
		h{k} = rand( size(x)-T+1);
		h{k} = h{k};
	end
    %% h is now a 1x88 cell, h{note} is a 5x100 array of random values.
    %% The 5 comes from the height of the CQ array minus the length of
    %% a template, plus 1. I guess this is space to allow for the
    %% 5-bins-per-semitone pitch shift.
end
if ~exist( 'z') || isempty( z)
    %% does happen, z was not provided
    z = cell(1, K);
	for k = 1:K
		z{k} = rand( size(x,2),1);
		z{k} = z{k};
	end
    %% z is a 1x88 cell, z{note} is a 100x1 array of random values.
end
if ~exist( 'u') || isempty( u)
    %% does happen, u was not provided
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
    %% u is a 10x88 cell, u{instrument,note} is a 100x1 double containing
    %% all 1s if note is in-range for instrument and all 0s otherwise
end

fh = cell(1, K); %% 1x88
fw = cell(R, K); %% 10x88
for k = 1:K
    fh{k} = ones(wc) + 1i*ones(wc);
    for r=1:R
        fw{r,k} = ones(wc) + 1i*ones(wc);
    end;
end;

%% now fh is a 1x88 cell, and fh{note} is a 553x199 array initialised
%% with all complex values 1 + 1i

%% fw is a 10x88 cell, and fw{instrument,note} is a 553x199 array
%% likewise


%% The MIREX paper describes the model as
%%
%% P(w,t) = P(t) sum[p,f,s]( P(w-f|s,p) P(f|p,t) P(s|p,t) P(p|t) )
%%
%% where 
%%   w = log-frequency bin index
%%   t = time
%%   s = instrument number
%%   p = MIDI pitch number
%%   f = frequency offset (pitch-adjustment convolution)
%% so
%%   P(w,t) = the input distribution (constant q spectrum)
%%   P(t) = overall energy by time
%%   P(w-f|s,p) = spectral template for instrument s, pitch p with offset f
%%   P(f|p,t) = frequency offset for pitch p, time t?
%%   P(s|p,t) = instrument contribution for pitch p, time t
%%   P(p|t) = pitch probability for time t
%% the outputs we want to produce are P(p|t) (the transcription matrix)
%% and P(s|p,t) (the instrument classification).
%%
%% In this program,
%%   x -> P(w,t), the input distribution
%%   w -> P(w|s,p), the templates
%%   h -> P(f|p,t), the pitch shift component
%%   z -> P(p|t), the pitch probabilities, the main return value
%%   u -> P(s|p,t), the source contribution, the secondary return value
%%
%% The paper gives the update rule for the expectation step as
%% 
%% P(p,f,s|w,t) =       P(w-f|s,p) P(f|p,t) P(s|p,t) P(p|t)
%%                --------------------------------------------------
%%                sum[p,f,s] ( P(w-f|s,p) P(f|p,t) P(s|p,t) P(p|t) )
%%
%% and the update equations for the maximization step as
%%
%% P(f|p,t) =  sum[w,s] ( P(p,f,s|w,t) P(w,t) )
%%  or h      ----------------------------------
%%            sum[f,w,s] ( P(p,f,s|w,t) P(w,t) )
%%
%% P(s|p,t) =  sum[w,f] ( P(p,f,s|w,t) P(w,t) )
%%  or u      ----------------------------------
%%            sum[s,w,f] ( P(p,f,s|w,t) P(w,t) )
%%
%% P(p|t) =  sum[w,f,s] ( P(p,f,s|w,t) P(w,t) )
%%  or z    ------------------------------------
%%          sum[p,w,f,s] ( P(p,f,s|w,t) P(w,t) )
%%
%% (there is also an update equation for x, or P(w|s,p) but we
%% don't want that as it's the input)



% Make commands for subsequent multidim operations and initialize fw

fnh = 'c(hc(1)-(T(1)+(1:size(h{k},1))-2),hc(2)-(T(2)+(1:size(h{k},2))-2))';
xai = 'xa(1:size(x,1),1:size(x,2))';
flz = 'xbar(end:-1:1,end:-1:1)';

for k = 1:K
    for r=1:R
        if( (pa(r,1) <= k &&  k <= pa(r,2)) )

	  %% fftn(X,siz) takes an N-dimensional FFT (same number of
	  %% dimensions as siz) of X, padding or truncating X
	  %% beforehand so that it is of size siz. Here w{r,k} is the
	  %% 545x1 template for instrument r and note k, and wc is
	  %% 553x199.

	  %% I believe this is equivalent to performing a 553-point
	  %% FFT of each column of the input (with w{r,k} in the first
	  %% 545 elements of the first column of that input) and then
	  %% a 199-point FFT of each row of the result.

	  %% The output is of course complex.

	  %% The purpose of this is to support convolution for pitch
	  %% shifting. w{r,k} are the templates, and fw{r,k} are ffts
	  %% of the templates which will be multiplied by fh, the
	  %% equivalent ffts of the pitch shift contributions, later

            fw{r,k} = fftn( w{r,k}, wc);
        end;
    end;
end;

% Iterate
for it = 1:iter
    
    %disp(['Iteration: ' num2str(it)]);
    
    % E-step
    xa = eps; %% tiny non-zero initialiser as we'll be dividing by this later
    for k = 16:73  %% overall note range found in instrument set
        fh{k} = fftn( h{k}, wc); %% this and the subsequent ifftn are for the pitch-shift convolution step I think
        for r=1:R  %% instruments
            if( (pa(r,1) <= k &&  k <= pa(r,2)) )
                xa1 = abs( real( ifftn( fw{r,k} .* fh{k})));                
                xa = xa + xa1(1:size(x,1),1:size(x,2)) .*repmat(z{k},1,size(x,1))'.*repmat(u{r,k},1,size(x,1))';

		%% so xa is the accumulation of the element-by-element
		%% product of: the pitch-shifted templates (xa1); the
		%% pitch probabilities (z); and the source
		%% contributions (u); across all instruments and notes

		%% note that xa1 is resized to match x, the input,
		%% which is possible because fw and fh were
		%% constructed at (just over) the same width and
		%% (almost) twice the height (553x199 if x is
		%% 549x100). 

		%% the other components are 100x1, i.e. one value per
		%% time step, so they are tiled up to 100x549 and then
		%% transposed for the multiplication

                clear xa1;
            end
        end
    end
    
    xbar = x ./ xa;
    xbar = eval( flz);
    fx = fftn( xbar, wc);
    
    
    % M-step
    for k = 16:73
        
        
      %% Throughout here, r is an instrument number (1..R) and k is a
      %% note number (1..K)

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
