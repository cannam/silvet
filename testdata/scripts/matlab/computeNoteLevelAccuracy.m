function [Pre,Rec,F,Acc,PreOff,RecOff,FOff,AccOff] = computeNoteLevelAccuracy(nmat1,nmat2)

% Compute note-level onset-only and onset-offset accuracy (Bay09)


% Initialize
if (isempty(nmat1)) Pre=0; Rec=0; F=0; Acc=0; return; end;

% Total number of transcribed notes
Ntot = size(nmat1,1);

% Number of reference notes
Nref = size(nmat2,1);

% Number of correctly transcribed notes, onset within a +/-50 ms range
Ncorr = 0;
NcorrOff = 0;
for j=1:size(nmat2,1)
    for i=1:size(nmat1,1)
        if( (nmat1(i,3) == nmat2(j,3)) && (abs(nmat2(j,1)-nmat1(i,1))<=0.05) )
            Ncorr = Ncorr+1;
            
            % If offset within a +/-50 ms range or within 20% of ground-truth note's duration
            if abs(nmat2(j,2) - nmat1(i,2)) <= max(0.05, 0.2 * (nmat2(j,2) - nmat2(j,1)))
               NcorrOff = NcorrOff +1;
            end;
            
            break; % In order to consider duplicates as false alarms
            
        end;
    end;
end;

% Number of onset-only P-R-F-Acc
Nfp = Ntot-Ncorr;
Nfn = Nref-Ncorr;
Rec = Ncorr/Nref;
Pre = Ncorr/Ntot;
F = 2*((Pre*Rec)/(Pre+Rec));
Acc= Ncorr/(Ncorr+Nfp+Nfn);

% Number of onset-offset P-R-F-Acc
NfpOff = Ntot-NcorrOff;
NfnOff = Nref-NcorrOff;
RecOff = NcorrOff/Nref;
PreOff = NcorrOff/Ntot;
FOff = 2*((PreOff*RecOff)/(PreOff+RecOff));
AccOff= NcorrOff/(NcorrOff+NfpOff+NfnOff);