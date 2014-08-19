function [wproj,wnorm] = compWprojW(W)
% [wproj,wnorm] = compWprojW(W)
%
% Computes the full time-convolution of each waveform with every other waveform (
% (This convolution is 'valid' in the rows, and 'full' in the columns)
%
% INPUT: 
%   W [ntime x nelectrodes x ncells ] - tensor of spike waveforms
% 
% OUTPUT: 
%   wproj [ 2*ntime-1 x ncells x ncells] - tensory of waveforms convolved w each other
%   wnorm [ ncells x 1 ] - squared L2 norm of each waveform
%
% jw pillow 8/18/2014


[nw,~,nc] = size(W);

wproj = zeros(2*nw-1,nc,nc); 
for jcell = 1:nc
    for icell = jcell:nc  % compute conv of w(icell) convolved with w(jcell)
        zz = W(:,:,jcell)*W(:,:,icell)';
        for j = 1:nw
            wproj(j:j+nw-1,icell,jcell) = wproj(j:j+nw-1,icell,jcell) + zz(:,nw-j+1);
        end
	% for all cell pairs (not including auto-correlation)
        if icell > jcell
            wproj(:,jcell,icell) = flipud(wproj(:,icell,jcell));
	end
    end
end
wnorm = diag(squeeze(wproj(nw,:,:)));  % dot prod of each waveform with itself
