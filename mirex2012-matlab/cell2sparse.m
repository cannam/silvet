function spCQT = cell2sparse(Xcq,octaves,bins,firstcenter,atomHOP,atomNr)
%Generates a sparse matrix containing the CQT coefficients (rasterized).
%
%The sparse matrix representation of the CQT coefficients contains all
%computed coefficients at the corresponding time-frequency location
%(similar to a spectrogram). For lower frequencies this means, that
%each coefficient is followed by zeros stemming from the fact, that the time
%resolution for lower frequencies decreases as the frequency resolution
%increases. Due to the design of the CQT kernel, however, the coefficients
%of different octaves are synchronised, meaning that for the second highest
%octave each coefficient is followed by one zero, for the next octave down
%two zeros are inserted, for the next octave four zeros are inserted and so
%on.
%
%INPUT:
%   Xcq         ... Cell array consisting of all coefficients for all octaves
%   octaves     ... Number of octaves processed
%   bins        ... Number of bins per octave
%   firstcenter ... Location of the leftmost atom-stack in the temporal
%                   kernel
%   atomHOP     ... Spacing of two consecutive atom stacks
%   atomNr      ... Number of atoms per bin within the kernel
%
%Christian Schörkhuber, Anssi Klapuri 2010-06

if 0
%% this version has big memory consumption but is very fast
    emptyHops = firstcenter/atomHOP;
    drop = emptyHops*2^(octaves-1)-emptyHops; %distance between first value in highest octave and first value in lowest octave
    spCQT = zeros(bins*octaves,size(Xcq{1},2)*atomNr-drop);

    for i=1:octaves
        drop = emptyHops*2^(octaves-i)-emptyHops; %first coefficients of all octaves have to be in synchrony
        X = Xcq{i}; 
        if  atomNr > 1 %more than one atom per bin --> reshape
            Xoct = zeros(bins,atomNr*size(X,2)-drop);
            for u=1:bins %reshape to continous windows for each bin (for the case of several wins per frame)
               octX_bin = X((u-1)*atomNr+1:u*atomNr,:);
               Xcont = reshape(octX_bin,1,size(octX_bin,1)*size(octX_bin,2));
               Xoct(u,:) = Xcont(1+drop:end);
            end
            X = Xoct;
        else
            X = X(:,1+drop:end);
        end
        binVec = bins*octaves-bins*i+1:bins*octaves-bins*(i-1);
        spCQT(binVec,1:2^(i-1):size(X,2)*2^(i-1)) = X;

    end
    spCQT = sparse(spCQT); %storing as sparse matrix at the end is the fastest way. Big memory consumption though!

else
%% this version uses less memory but is noticable slower
    emptyHops = firstcenter/atomHOP;
    drops = emptyHops*2.^(octaves-(1:octaves))-emptyHops;
    len = max(((atomNr*cellfun('size',Xcq,2)-drops).*2.^(0:octaves-1))); %number of columns of output matrix
    spCQT = [];

    for i=octaves:-1:1
        drop = emptyHops*2^(octaves-i)-emptyHops; %first coefficients of all octaves have to be in synchrony
        X = Xcq{i}; 
        if  atomNr > 1 %more than one atom per bin --> reshape
            Xoct = zeros(bins,atomNr*size(X,2)-drop);
            for u=1:bins %reshape to continous windows for each bin (for the case of several wins per frame)
               octX_bin = X((u-1)*atomNr+1:u*atomNr,:);
               Xcont = reshape(octX_bin,1,size(octX_bin,1)*size(octX_bin,2));
               Xoct(u,:) = Xcont(1+drop:end);
            end
            X = Xoct;
        else
            X = X(:,1+drop:end);
        end
        X = upsample(X.',2^(i-1)).';
        X = [X zeros(bins,len-size(X,2))];
        spCQT = sparse([spCQT; X]);  
    end

end
