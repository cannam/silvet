function  [Pre,Rec,F]  = batchProcessingEvaluate(folder)

% Evaluate transcription output, using onset-only note-based metrics
% e.g. [Pre,Rec,F]  = batchProcessingEvaluate('TRIOS-mirex2012-matlab');


fileList = dir(folder);
fileCount = 0;

for i=3:length(fileList) 
    
    if(isdir([folder '/' fileList(i).name])) 
        
        fileCount = fileCount + 1;

        
        % Load ground truth
        [pianoRollGT,nmatGT] = convertMIDIToPianoRoll([fileList(i).name '.mid'],10,1.0);
        
        
        % Load transcripton nmat
        nmat = load([folder '/' fileList(i).name '/' 'mix.lab']);
        
        
        % Convert 3rd nmat column to MIDI scale        
        nmat(:,3) = round(12.*log2(nmat(:,3)./27.5) + 1);
        
        
        % Compute onset-based note-level accuracy
        [Pre(fileCount),Rec(fileCount),F(fileCount)] = computeNoteLevelAccuracy(nmat,nmatGT);
        
    end; 
end;