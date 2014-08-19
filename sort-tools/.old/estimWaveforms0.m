function What = estimWaveforms0(yy, twin, Xsp, nw, CriticalSize);
%  What = estimWaveforms0(yy, twin, Xsp, nw, CriticalSize);
%
%  Computes estimate of spike waveform for cells from spike trains and
%  electrode data, using (correct) least-squares regression
%
%  Input:  
%   yy = electrode data;
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
Xxc = zeros(2*nw-1,nc,nc);
XY = zeros(nw*nc,ne);

% Decide how many chunks to use
i1 = istrt;
maxchun = round(CriticalSize/ne);
i2 = min(i1+maxchun, iend);
fprintf(1, 'Computing waveforms in %d steps\n', ceil(slen/maxchun));

% ---------------------------------------------------
% 1. Compute correlation of X with itself, and X with Y
while i1 < iend
    fprintf('step num: %d, starting at %d\n', round(i1/maxchun), i1);
    ywin = [i1-ipre-1, i2-ipre+nw-1]+twin(1);
    Ydat = yy(ywin(1)+1:ywin(2),:);
    for jcell = 1:nc
        iisp = find(Xsp(i1:i2,jcell));
        xx0 = sparse(2*nw-1,nc);
        yy0 = zeros(nw,ne);
        for jsp = 1:length(iisp);
            xx0 = xx0+Xsp(i1+iisp(jsp)-nw:i1+iisp(jsp)+nw-2,:);
            yy0 = yy0+Ydat(iisp(jsp):iisp(jsp)+nw-1,:);
        end
        Xxc(:,:,jcell) = Xxc(:,:,jcell)+xx0;
        XY((jcell-1)*nw+1:jcell*nw,:) = XY((jcell-1)*nw+1:jcell*nw,:) + yy0;
    end
    i1 = i2+1;
    i2 = min(i1+maxchun,iend);
end

% ---------------------------------------------------
% 2. Insert cross-corrs into giant covariance matrix for X'X.
Xcov = zeros(nw*nc);
ii1 = [nw:-1:1];
ii2 = nw:2*nw-1;
for j = 1:nc
    for i = j:nc;
        ncol = (j-1)*nc+i;
        if i == j
            Xcov((i-1)*nw+1:i*nw,(j-1)*nw+1:j*nw) = ...
                toeplitz(Xxc(ii2,j,i));
        else
            aa = toeplitz(Xxc(ii1,j,i),Xxc(ii2,j,i));
            Xcov((i-1)*nw+1:i*nw,(j-1)*nw+1:j*nw) = aa';
            Xcov((j-1)*nw+1:j*nw,(i-1)*nw+1:i*nw) = aa;
        end
    end
end

% ----------------------------------------------------
% 3. Compute least-squares solution; reshape matrix
ww = Xcov\XY;
for i = 1:nc
    What(:,:,i) = ww((i-1)*nw+1:i*nw,:);
end
