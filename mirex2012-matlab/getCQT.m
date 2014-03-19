function intCQT = getCQT(Xcqt,fSlice,tSlice,iFlag)
%outCQ = getCQT(Xcqt,fSlice,tSlice,iFlag) computes a rasterized representation of 
%the amplitudes of the calculated CQT coefficients for the frequency bins definded in vector fSlice and the 
%points in time (time frames) defined in vector tSlice using the interpolation method defined in iFlag. 
%Valid values for iFlag are:
%
%'linear'  ... linear interpolation (default)
%'spline'  ... spline interpolation
%'nearest' ... nearest neighbor interpolation
%'cubic'   ... piecewise cubic interpolation
%
%If the entire CQT representation should be rasterized, set fSlice and
%tSlice to 'all'.
%The input parameter Xcqt is the structure gained using cqt(...).
%The output parameter 'intCQT' is the same size as Xcqt.spCQT but is no
%longer sparse since the zeros between two coefficients are replaced by
%the interpolated values. The coefficients stored in 'intCQT' are now
%real-valued since only the absolute values of the coefficients are
%interpolated. If a spectrogram-like (rasterized) version of the CQT
%coefficients including phase information is required, use the function
%cqtPerfectRast() (see documentation for further information)
%
%Christian Schörkhuber, Anssi Klapuri 2010-06


if ischar(fSlice), fSlice = 1:(Xcqt.bins*Xcqt.octaveNr); end;
if ischar(tSlice)
    lastEnt = find(Xcqt.spCQT(1,:),1,'last');
    tSlice = 1:lastEnt;
end
if nargin < 4, iFlag = 'linear'; end;

intCQT = zeros(length(fSlice),length(tSlice));
bins = Xcqt.bins;
spCQT = Xcqt.spCQT;
octaveNr = Xcqt.octaveNr;
spCQT = spCQT.';

for k=1:length(fSlice)
   oct = octaveNr-floor((fSlice(k)-0.1)/bins);
   stepVec = 1:2^(oct-1):size(spCQT,1);
   Xbin = full(spCQT(stepVec,fSlice(k)));
   intCQT(k,:) = interp1(stepVec,abs(Xbin),tSlice,iFlag);
end



