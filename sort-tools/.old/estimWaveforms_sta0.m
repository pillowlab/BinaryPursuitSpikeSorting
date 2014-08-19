function What = estimWaveforms_sta0(yy, twin, Xsp, nw, CriticalSize);
%  What = estimWaveforms_sta0(yy, twin, Xsp, nw, CriticalSize);
%
%  Computes estimate of spike waveform for cells from spike trains and
%  electrode data, using just the AVERAGE waveform associated w/ a spike
%
%  Input:  
%   yy = electrode data
%   twin = time window;
%   Xsp = spikes;  coloumn vector of spike count in each time bin
%   nw  = number of time samples to include in the spike-triggered ensemble
%   CriticalSize - maximum number of floats to store while computing cov
%                  (smaller = slower but smaller memory requirement)
%
%  Output:  What = estimated waveform for each cell


% maximum size of chunks -- decrease this if getting "out of memory" errors
if nargin < 5
    CriticalSize = 1e7;
end

[slen,ne] = size(yy);  % electrode data
[slen,nc] = size(Xsp); % spike data

Xsp = [zeros(nw/2-1,nc); Xsp; zeros(nw/2,nc)]; % pad edges with zeros
ipre = nw-1;
istrt = nw;
iend = slen;

% --------------------
What = zeros(nw,ne,nc);

% Decide how many chunks to use
i1 = istrt;
maxchun = round(CriticalSize/ne);
i2 = min(i1+maxchun, iend);
fprintf(1, 'Computing waveforms in %d steps\n', ceil(slen/maxchun));

% ---------------------------------------------------
% 1. Compute correlation of X with Y
while i1 < iend
    fprintf('step num: %d, starting at %d\n', round(i1/maxchun), i1);
    ywin = [i1-ipre-1, i2-ipre+nw-1]+twin(1);
    Ydat = yy(ywin(1)+1:ywin(2),:);
    for jcell = 1:nc
        iisp = find(Xsp(i1:i2,jcell)); 
        nsp = length(iisp);
        yy0 = zeros(nw,ne);
        for jsp = 1:length(iisp);
            yy0 = yy0+Ydat(iisp(jsp):iisp(jsp)+nw-1,:);
        end
        What(:,:,jcell) = yy0/nsp;
    end
    i1 = i2+1;
    i2 = min(i1+maxchun,iend);
end

