function [What,wwsigs]=estimWaveforms(X, Y, nw)
%  [What,wwsigs]=estimWaveforms(X, Y, nw)
%
%  Computes estimate of spike waveform for cells from spike trains and
%  electrode data, using (correct) least-squares regression
%
%  Input:  
%  ------
%   X [nsamps x ncells] - each column holds spike train of single neuron
%   Y [nsamps x nelec] - raw electrode data
%   nw [1 x 1] - number of time bins in the spike waveform
%
%  Output:  
%  -------
%   What [nw x ne x ncells] - estimated waveforms for each cell
%   sigs [ncells x 1] - posterior stdev of each neuron's waveform coefficients
%
% jw pillow 8/18/2014


[nt,nc] = size(X); % number of time bins and number of cells
ne = size(Y,2); % number of electrodes
nw2 = nw/2;

% Compute blocks for covariance matrix XX and cross-covariance XY
XXblocks = zeros(nc*nw,nc);
XY = zeros(nc*nw,ne);
for jj = 1:nw
    inds = ((jj-1)*nc+1):(jj*nc);
    XXblocks(inds,:) = X(1:end-jj+1,:)'*X(jj:end,:); % spike train covariance
    XY(inds,:) = X(max(1,nw2-jj+2):min(nt,nt-jj+nw2+1),:)'*... 
        Y(max(1,jj-nw2):min(nt,nt+jj-nw2-1),:); % cross-covariance 
end

% Insert blocks into covariance matrix XX
XX = zeros(nc*nw,nc*nw);
for jj = 1:nw
    inds1 = ((jj-1)*nc+1):(nc*nw);
    inds2 = ((jj-1)*nc+1):(jj*nc);
    XX(inds1,inds2) = XXblocks(1:(nw-jj+1)*nc,:);  % below diagonal blocks
    XX(inds2,inds1) = XXblocks(1:(nw-jj+1)*nc,:)'; % above diagonal blocks
end
What = XX\XY; % do regression
What = permute(reshape(What,[nc,nw,ne]),[2,3,1]);  % reshape into tensor

% 4. If desired, compute posterior variance for each waveform (function of # spikes)
if nargout > 1
    wwsigs = sqrt(1./diag(XX(1:nc,1:nc)));
end
 % Note: the "correct" formula should be diag(inv(Xcov)), but this is
 % close, ignoring edge effects, and much faster;