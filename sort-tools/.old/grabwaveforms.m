function ysamps = grabwaveforms(Ydat,sps,nw);
% ysamps = grabwaveforms(Ydat,sps,nw);
%
% Grabs the raw spike waveform associated with each spike in sps;

wwid = nw/2;
spinds = find(sps);
spinds = spinds((spinds>wwid) & (spinds<=size(Ydat,1)-wwid));
nsp = length(spinds);
ne = size(Ydat,2);
ysamps =  zeros(nsp,ne*nw);

for j = 1:length(spinds)
    ysmp = Ydat(spinds(j)-wwid:spinds(j)+wwid-1,:);
    ysamps(j,:) = ysmp(:)';
end
    

