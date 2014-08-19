function Vpred  = compVpredictionSprse(xsp, W)
% Vpred  = compVpredictionSprse(xsp, ww)
%
% Computes sparse binary xsp convolved with ww spike waveforms ww
% 
% INPUT:  
%   xsp [nt x nneur] - sparse binary matrix, each column is a spike train
%   ww [nw x nelec x nneur] - tensor of spike waveforms
%
% OUTPUT:
%   Vpred [nt x nelec] - convolution of xsp with ww
%
% jw pillow 8/18/2014


nc = size(xsp,2);  % number of cells
wwid = size(W,1)/2;  % 1/2 length of spike waveform
iirel = (0:wwid*2-1)'; % relative time indices
slen = size(xsp,1);  % number of time samples
  
Vpred = zeros(slen+wwid*2,size(W,2));  % allocate memory (with padding at beginning and end)

for j = 1:nc  % loop over neurons
    isp = find(xsp(:,j));  % find the spike times
    for i = 1:length(isp);
        ii = isp(i)+iirel;
        Vpred(ii,:) = Vpred(ii,:)+W(:,:,j);  
    end
end
Vpred = Vpred(wwid+1:end-wwid,:);  % remove padding at the end

