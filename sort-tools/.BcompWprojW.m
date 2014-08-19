function [wproj,wnorm] = compWprojW(WW,eeinter);
% [wproj,wnorm] = compWprojW(WW,eeintersect);
%
% Computes the full time-convolution of each waveform with every other waveform

[nw,ne,nc] = size(WW);

wproj = zeros(2*nw-1,nc,nc);
for jcell = 1:nc
    for icell = jcell:nc
        zz = WW(:,eeinter{jcell}{icell},jcell)*WW(:,eeinter{jcell}{icell},icell)';
        for j = 1:nw
            wproj(j:j+nw-1,icell,jcell) = wproj(j:j+nw-1,icell,jcell) + zz(:,nw-j+1);
        end
        if icell > jcell
            wproj(:,jcell,icell) = flipud(wproj(:,icell,jcell));
        end
    end
end
wnorm = diag(squeeze(wproj(nw,:,:)));  % dot prod of each waveform with itself


% % -- Older code (need wwe instead of eeintersect) ----
% for jcell = 1:10
%    wproj(:,:,jcell) = fullfilt(WW(:,wwe{jcell},jcell), WW(:,wwe{jcell},:));
% end
