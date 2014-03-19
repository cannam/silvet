function Xcqt = cqt(x,fmin,fmax,bins,fs,varargin) 
%Xcqt = cqt(x,fmin,fmax,bins,fs,varargin) calculates the constant-Q transform of the input signal x.
%
%INPUT:
%   fmin      ... lowest frequency of interest
%   fmax      ... highest frequency of interest
%   bins      ... frequency bins per octave
%   fs        ... sampling rate
%
%   optional input parameters (parameter name/value pairs):
%
%   'atomHopFactor'   ... overlap of temporal atoms in percent. Default: 0.25.
%    
%   'q'       ... the maximum value for optimal reconstruction is q=1.
%                 For values smaller than 1 the bandwidths of the spectral
%                 atoms (filter) are increased retaining their center
%                 frequencies (frequency 'smearing', frequency domain redundancy 
%                 increases, time resolutin improves). Default: 1.
%   'thresh'  ... all values in the cqt kernel smaller than tresh are
%                 rounded to zero. A high value for thresh yields a
%                 very sparse kernel (fast) but introduces a bigger error. 
%                 The default value is chosen so that the error due to rounding is negligible.
%   'kernel'  ... if the cqt kernel structure has been precomputed
%                 (using function 'genCQTkernel'), the computation of the kernel
%                 will be by-passed below).
%   'win'     ... defines which window will be used for the CQT. Valid
%                 values are: 'blackman','hann' and 'blackmanharris'. To
%                 use the square root of each window use the prefix 'sqrt_'
%                 (i.e. 'sqrt_blackman'). Default: 'sqrt_blackmanharris'
%   'coeffB',
%   'coeffA'  ... Filter coefficients for the anti-aliasing filter, where
%                 'coeffB' is the numerator and 'coeffA' is the
%                 denominator (listed in descending powers of z). 
%                                                  
%OUTPUT:
%   Xcqt      ... struct that comprises various fields: 
%              spCQT: CQT coefficients in the form of a sparse matrix 
%                    (rasterized, not interpolated)
%              fKernel: spectral Kernel 
%              fmin: frequency of the lowest bin
%              fmax: frequency of the hiqhest bin
%              octaveNr: number of octaves processed
%              bins: number of bins per octave
%              intParams: structure containing additional parameters for the inverse transform   
%
%Christian Schörkhuber, Anssi Klapuri 2010-06

%% input checking
if size(x,2) > 1 && size(x,1) > 1, error('cqt requires one-dimensional input!'); end;
if size(x,2) > 1, x = x(:); end; %column vector

%% input parameters
q = 1; %default value
atomHopFactor = 0.25; %default value
thresh = 0.0005; %default value
winFlag = 'sqrt_blackmanharris';

for ain = 1:1:length(varargin)
    if strcmp(varargin{ain},'q'), q = varargin{ain+1}; end;
    if strcmp(varargin{ain},'atomHopFactor'), atomHopFactor = varargin{ain+1}; end;
    if strcmp(varargin{ain},'thresh'), thresh = varargin{ain+1}; end;
    if strcmp(varargin{ain},'kernel'), cqtKernel = varargin{ain+1}; end;
    if strcmp(varargin{ain},'win'), winFlag = varargin{ain+1}; end;
    if strcmp(varargin{ain},'coeffB'), B = varargin{ain+1}; end;
    if strcmp(varargin{ain},'coeffA'), A = varargin{ain+1}; end;    
end

%% define
octaveNr = ceil(log2(fmax/fmin));
xlen_init = length(x);

%% design lowpass filter
if ~exist('B','var') || ~exist('A','var')
    LPorder = 6; %order of the anti-aliasing filter
    cutoff = 0.5;
    [B A] = butter(LPorder,cutoff,'low'); %design f_nyquist/2-lowpass filter
end

%% design kernel for one octave 
if ~exist('cqtKernel','var')
    cqtKernel = genCQTkernel(fmax, bins,fs,'q',q,'atomHopFactor',atomHopFactor,'thresh',thresh,'win',winFlag);
end

%% calculate CQT
cellCQT = cell(1,octaveNr);
maxBlock = cqtKernel.fftLEN * 2^(octaveNr-1); %largest FFT Block (virtual)
suffixZeros = maxBlock; 
prefixZeros = maxBlock;
x = [zeros(prefixZeros,1); x; zeros(suffixZeros,1)]; %zeropadding
OVRLP = cqtKernel.fftLEN - cqtKernel.fftHOP;
K = cqtKernel.fKernel'; %conjugate spectral kernel for cqt transformation

for i = 1:octaveNr
    xx = buffer(x,cqtKernel.fftLEN, OVRLP,'nodelay'); %generating FFT blocks
    XX = fft(xx); %applying fft to each column (each FFT frame)
    cellCQT{i} = K*XX; %calculating cqt coefficients for all FFT frames for this octave

    if i~=octaveNr
        x = filtfilt(B,A,x); %anti aliasing filter
        x = x(1:2:end); %drop samplerate by 2
    end
end

%% map to sparse matrix representation
spCQT = cell2sparse(cellCQT,octaveNr,bins,cqtKernel.firstcenter,cqtKernel.atomHOP,cqtKernel.atomNr);

%% return
intParam = struct('sufZeros',suffixZeros,'preZeros',prefixZeros,'xlen_init',xlen_init,'fftLEN',cqtKernel.fftLEN,'fftHOP',cqtKernel.fftHOP,...
    'q',q,'filtCoeffA',A,'filtCoeffB',B,'firstcenter',cqtKernel.firstcenter,'atomHOP',cqtKernel.atomHOP,...
    'atomNr',cqtKernel.atomNr,'Nk_max',cqtKernel.Nk_max,'Q',cqtKernel.Q,'rast',0);

Xcqt = struct('spCQT',spCQT,'fKernel',cqtKernel.fKernel,'fmax',fmax,'fmin',fmin*2^(1/bins),'octaveNr',octaveNr,'bins',cqtKernel.bins,'intParams',intParam);


