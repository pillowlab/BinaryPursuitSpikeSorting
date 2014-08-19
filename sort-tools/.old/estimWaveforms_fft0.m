function What = estimWaveforms_fft0(yy, twin, xsp, nw, CriticalSize);
%  What = estimWaveforms_fft0(yy, twin, xsp, nw, CriticalSize);
%
%  Computes estimate of spike waveform for cells whose spike trains
%
%  Input:  
%   yy = electrode data
%   twin = time window;
%   xsp = spikes;  coloumn vector of spike count in each time bin
%   n  = number of time samples to include in the spike-triggered ensemble
%   CriticalSize - maximum number of floats to store while computing cov
%                    (smaller = slower but smaller memory requirement) 
%
%  Output:  What = estimated waveform for each cell


% maximum size of chunks -- decrease this if getting "out of memory" errors
if nargin < 5
    CriticalSize = 1e8;
end
[slen,ne] = size(yy);
[slen,nc] = size(xsp);

% -----------------------------------------------
% 1. Compute cross-correlations between xsp and itself
xxc = zeros(nw*2-1,nc,nc);
xxhat = fft(full(xsp));
ii0 = 1:nw;
iikp = [slen-nw+2:slen, 1:nw];

for j = 1:nc
    for i = j:nc
        if i == j
            q = ifft(abs(xxhat(:,j)).^2);
            xxc(ii0,j,i) = q(ii0);
        else
            q = ifft(conj(xxhat(:,j)).*xxhat(:,i));
            xxc(:,j,i) = q(iikp);
            %xxc(:,i,j) = flipud(xxc(:,j,i));
      end
    end
end

% ---------------------------------------------------
% 2. Insert cross-corrs into giant covariance matrix
% Insert elements of xxc into matrix for YX and 
ii1 = [nw:-1:1];
ii2 = nw:2*nw-1;
for j = 1:nc
    for i = j:nc;
        ncol = (j-1)*nc+i;
        if i == j
            xx((i-1)*nw+1:i*nw,(j-1)*nw+1:j*nw) = ...
                toeplitz(xxc(ii0,j,i));
        else
            aa = toeplitz(xxc(ii1,j,i),xxc(ii2,j,i));
            xx((i-1)*nw+1:i*nw,(j-1)*nw+1:j*nw) = aa;
            xx((j-1)*nw+1:j*nw,(i-1)*nw+1:i*nw) = aa';
        end
    end
end

% ---------------------------------------------------
% 3. Compute cross-corr of x with y, insert into 'xy'
xy = zeros(nw*nc,ne);
iikp2 = [slen-nw/2+1:slen, 1:nw/2];
for j = 1:ne
    yh = fft(yy(:,j));
    for i = 1:nc
        q = ifft(conj(xxhat(:,i)).*yh);
        xy((i-1)*nw+1:i*nw,j) = q(iikp2);
    end
end

ww = xx\xy;
for i = 1:nc
    What(:,:,i) = ww((i-1)*nw+1:i*nw,:);
end

