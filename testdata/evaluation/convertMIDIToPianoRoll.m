function [pianoRoll,newnmat] = convertMIDIToPianoRoll(filename,timeResolution,dur)

% Time resolution is in msec
% eg. pianoRoll = convertMIDIToPianoRoll('bach_847MINp_align.mid',10);

% Read MIDI file
nmat = readmidi(filename);


% Gather MIDI information
[n1 n2] = size(nmat);
lenthInSec = nmat(n1,6) + nmat(n1,7);
pianoRoll = zeros(88,round(lenthInSec*(1000/timeResolution)));


% Fill piano roll
for i=1:n1
    pianoRoll(round(nmat(i,4)-20),round(nmat(i,6)*(1000/timeResolution))+1:round(nmat(i,6)*(1000/timeResolution)+dur*nmat(i,7)*(1000/timeResolution))+1) = 1;
end;

pianoRoll = pianoRoll';

% Plot piano roll
%figure; imagesc(imrotate(pianoRoll,0)); axis xy
%colormap('gray'); xlabel('time frame'); ylabel('Pitch');


% Convert to non-MIDI nmat
newnmat(:,1) = nmat(:,6);
newnmat(:,2) = nmat(:,6)+nmat(:,7);
newnmat(:,3) = nmat(:,4)-20;