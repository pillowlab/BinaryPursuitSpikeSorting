function [Ywht, tfilts,xfilts] = compWhiteningFilts(X,Y,W,nxc_t,nxc_x)
% [tfilts,xfilts] = compWhiteningFilts(X,Y,W,nxc_t,nxc_x)
%
% Compute filters that whitens y residuals
%
% INPUTS:
% ------
%  X [nsamps x ncells] - each column holds spike train of single neuronjj
%  Y [nsamps x nelec] - raw electrode data
%  W [nw x nelec x ncells] - tensor of spike waveforms
%  nxc_t [1x1]
%  nxc_x [1x1]
%
% OUTPUTS:
% --------
%  Ywht [nsamps x nelec] - whitened electrode data
%  tfilts [nxc_t x ncells] - temporal whitening filters
%  xfilts [nxc_x x nelec x nelec] spatial whitening filters
% 
%
%  NOTES:
% "Full" covariance too large to invert, so instead we proceed by:
%  1. Whitening each electrode in time (using nxc_t-tap filter)
%  2. Whiten across electrodes (using nxc_x timebins per electrode)
%
%  Algorithm:  
%    1. Compute the electrode residuals (raw electrodes minus spike
%       train convolved with waveforms) 
%    2. Compute auto-correlations for each electrode
%    3. Solve for temporal-whitening filter for each electrode
%    4. Temporally whiten electrode residuals
%    5. Compute spatially-whitening filters

slen = size(X,1); % spike train data (sparse)
ne = size(W,2); % number of electrodes

% === 1. Compute temporal whitening filts =================
ypred = compVpredictionSprse(X,W); % predicted electrode data based on spikes and W
yresid = Y-ypred;  % residuals 

% Method 1 (cost indep of nxc_t, and faster than xcorr)
tfilts = zeros(nxc_t,ne); % temporal whitening filters
for j = 1:ne
    % Compute autocovariance
    xc = circxcorr(yresid(:,j),nxc_t-1,'none');
    yxc = flipud(xc(1:nxc_t))/slen; % autocovariance

    % Compute whitening filt
    M = sqrtm(inv(toeplitz(yxc))); % whitening matrix
    tfilts(:,j) = M(:,nxc_t/2);
        
end
% % Check that it worked:
% plot(xcorr(conv2(yresid(:,j),tfilts(:,j),'valid'),nxc_t)/slen);


%  === 2. Compute spatial whitening filters ====================
% Will be faster when we used only neighboring electrodes)

% Compute temporally whiten residuals
yresid_wht = zeros(slen,ne);  % whitened residuals
for j = 1:ne; 
    yresid_wht(:,j) = conv2(yresid(:,j), tfilts(:,j), 'same');
end

% Compute spatial cross-correlation(s)
xxc = zeros(ne,ne,nxc_x);
xxc(:,:,1) = yresid_wht'*yresid_wht;
for j = 2:nxc_x
    xxc(:,:,j) = circshift(yresid_wht,j-1)'*yresid_wht;
end
xxc = xxc/slen;

% Insert into big covariance matrix 
M = zeros(ne*nxc_x);
for j = 1:ne
    for i = j:ne
        if i == j
            jj = (j-1)*nxc_x+1:j*nxc_x;
            M(jj,jj) = toeplitz(squeeze(xxc(j,j,:)));
        else
            jj = (j-1)*nxc_x+1:j*nxc_x;
            ii = (i-1)*nxc_x+1:i*nxc_x;
            M(jj,ii) = toeplitz(squeeze(xxc(j,i,:)),squeeze(xxc(i,j,:)));
            M(ii,jj) = M(jj,ii)';
        end
    end
end
% Compute filters
Q = sqrtm(inv(M));
xfilts = zeros(nxc_x,ne,ne);
for j = 1:ne
    xfilts(:,:,j) = reshape(Q(:,(j-1)*nxc_x+ceil(nxc_x/2)),[],ne);
end

% % ====================================
% % OPTIONAL: check that it worked
% for j = 1:ne
%     ywht(:,j) = conv2(yresid_wht,fliplr(xfilts(:,:,j)),'valid');
% end
% subplot(121);
% imagesc(cov(ywht));
% subplot(122);
% plot(xcorr(ywht(:,1:2),5));
% % ====================================

% === 3. Now whiten the raw data ==========================
Ywht = zeros(slen,ne);
for ii = 1:ne  % Temporally whiten
    Ywht(:,ii) = conv2(Y(:,ii),tfilts(:,ii),'same');
end;
Ywht = samefilt(Ywht,xfilts,'conv');  % spatially whiten
