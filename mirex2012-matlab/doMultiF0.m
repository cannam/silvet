function  []  = doMultiF0(inputFile,outputFile)


% Transcribe file
fprintf('%s',['Preprocessing............']);
[ph pz sumY] = transcriptionMultipleTemplates(inputFile,12,1.1,1.3);
fprintf('\n');
fprintf('%s',['Postprocessing...........']);
pianoRoll = repmat(sumY,88,1).*pz(1:88,:);
pianoRoll = pianoRoll';
for j=[1:15 74:88] pianoRoll(:,j)=0; end;
pianoRoll = medfilt1(pianoRoll,3);


% Max polyphony 5
[B,IX] =sort(pianoRoll,2,'descend');  
tempPianoRoll = zeros(size(pianoRoll,1),88);
for j=1:size(pianoRoll,1) for k=1:5 tempPianoRoll(j,IX(j,k)) = B(j,k); end; end;
pianoRoll = tempPianoRoll;


% Expand piano-roll and perform thresholding
expandedPianoRoll = zeros(4*size(pianoRoll,1),88);
for j=1:4*size(pianoRoll,1)
    expandedPianoRoll(j,:) = pianoRoll(floor((j-1)/4)+1,:);    
end;
finalPianoRoll = (expandedPianoRoll>4.8)';


% Create output and perform minimum duration pruning
auxPianoRoll = diff([zeros(1,88); finalPianoRoll'; zeros(1,88);],1); k=0;
for i=1:88
    onsets = find(auxPianoRoll(:,i)==1);
    offsets = find(auxPianoRoll(:,i)==-1);
    for j=1:length(onsets)           
        if((offsets(j)/100-0.01) - (onsets(j)/100) > 0.05)
            k=k+1;
            nmat(k,1) = onsets(j)/100;       
            nmat(k,2) = offsets(j)/100-0.01;    
            nmat(k,3) = 27.5*2.^((( i-1)*10 )/120);
        end;
    end;
end;
nmat = sortrows(nmat,1);


% Print output
fid=fopen(outputFile,'w');
for i=1:size(nmat,1)
    fprintf(fid,'%.2f\t%.2f\t%.2f\n',nmat(i,1),nmat(i,2),nmat(i,3));
end;
fclose(fid);
fprintf('%s','done');
fprintf('\n');


exit;
