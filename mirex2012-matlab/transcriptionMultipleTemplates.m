function [ph pz sumY] = transcriptionMultipleTemplates(filename,iter,sz,su)


% Load note templates
load('noteTemplatesBassoon');      W(:,:,1) = noteTemplatesBassoon;
load('noteTemplatesCello');        W(:,:,2) = noteTemplatesCello;
load('noteTemplatesClarinet');     W(:,:,3) = noteTemplatesClarinet;
load('noteTemplatesFlute');        W(:,:,4) = noteTemplatesFlute;
load('noteTemplatesGuitar');       W(:,:,5) = noteTemplatesGuitar;
load('noteTemplatesHorn');         W(:,:,6) = noteTemplatesHorn;
load('noteTemplatesOboe');         W(:,:,7) = noteTemplatesOboe;
load('noteTemplatesTenorSax');     W(:,:,8) = noteTemplatesTenorSax;
load('noteTemplatesViolin');       W(:,:,9) = noteTemplatesViolin;
load('noteTemplatesSptkBGCl');     W(:,:,10) = noteTemplatesSptkBGCl;


%pitchActivity = [14 16 30 40 20 21 38 24 35 1; 52 61 69 76 56 57 71 55 80 88]';
pitchActivity = [16 16 30 40 20 21 38 24 35 16; 52 61 69 73 56 57 71 55 73 73]';


W = permute(W,[2 1 3]);
W0 = squeeze(num2cell(W,1))';
clear('noteTemplatesBassoon','noteTemplatesCello','noteTemplatesClarinet','noteTemplatesFlute','noteTemplatesGuitar',...
    'noteTemplatesHorn','noteTemplatesOboe','noteTemplatesTenorSax','noteTemplatesViolin','noteTemplatesSptkBGCl','W');


% Compute CQT
[intCQT] = computeCQT(filename);
X = intCQT(:,round(1:7.1128:size(intCQT,2)))';
noiseLevel1 = medfilt1(X',40);
noiseLevel2 = medfilt1(min(X',noiseLevel1),40);
X = max(X-noiseLevel2',0);
Y = X(1:4:size(X,1),:);  % 40ms step
sumY = sum(Y');
clear('intCQT','X','noiseLevel1','noiseLevel2');

fprintf('%s','done');
fprintf('\n');
fprintf('%s',['Estimating F0s...........']);

% For each 2sec segment, perform SIPLCA with fixed W0
ph = zeros(440,size(Y,1));
pz = zeros(88,size(Y,1));

for j=1:floor(size(Y,1)/100)

    x=[zeros(2,100); Y(1+(j-1)*100:j*100,:)'; zeros(2,100)];
    [w,h,z,u,xa] = cplcaMT( x, 88, [545 1], 10, W0, [], [], [], iter, 1, 1, sz, su, 0, 1, 1, 1, pitchActivity);
   
    H=[]; for i=1:88 H=[H; h{i}]; end;
    ph(:,1+(j-1)*100:j*100) = H;
    Z=[]; for i=1:88 Z=[Z z{i}]; end;
    pz(:,1+(j-1)*100:j*100) = Z';   
    perc = 100*(j/(floor(size(Y,1)/100)+1));
    fprintf('\n');
    fprintf('%.2f%% complete',perc);
end;

len=size(Y,1)-j*100; % Final segment

if (len >0)
x=[zeros(2,len); Y(1+j*100:end,:)'; zeros(2,len)];
[w,h,z,u,xa] = cplcaMT( x, 88, [545 1], 10, W0, [], [], [], iter, 1, 1, sz, su, 0, 1, 1, 1, pitchActivity);
fprintf('\n');
fprintf('100%% complete');

H=[]; for i=1:88 H=[H; h{i}]; end;
ph(:,1+j*100:end) = H;
Z=[]; for i=1:88 Z=[Z z{i}]; end;
pz(:,1+j*100:end) = Z'; 
end;