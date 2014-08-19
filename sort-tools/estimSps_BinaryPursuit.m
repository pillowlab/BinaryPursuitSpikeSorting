function [Xhat,nrefviols] = estimSps_BinaryPursuit(loadYfun,wfile,nperblck,X0,pspike,nperchnk,minISI)
% [Xhat,nrefviols] = estimSps_BinaryPursuit(loadYfun,wfile,nsampsPerBlock,X0,pspike,blksize,minISI)
%
%  Estimates spike train over a large chunk of data by calling
%  runBinaryPursuit.m on many, smaller chunks
%
%  INPUTS
%     loadYfun - function for loading electrode data (takes single argument 'ywin')
%     wfile - filename for spike waveforms (with '%d' for block #)
%     nperblck - number of samples per block
%     X0 [nsamps x ncells] - initial guess at spike trains (sparse)
%     pspike - prior probability of a spike in a bin (pos scalar < 1)
%     nperchnk - numer of samples in a "chunk" to process at once for BP
%     minISI - minimum ISI in # samples
%
% OUTPUTS
%   Xhat = estimated binary spike train
%   nrefviols - number of spikes removed (post hoc) due to refractory violation
%
% jw pillow 8/18/2014


verbose = 10;  % report progress mod this value

[slen,nc] = size(X0); % total size of spike train data matrix
nblocks = slen/nperblck; % number of blocks (each with distinct waveform estimate)
nchnks = slen/nperchnk; % number of chunks to use for processing data
nchnksperblock = nperblck/nperchnk;

% Load 1st-block spike waveforms
wnum = 1; 
load(sprintf(wfile, wnum),'W');
nw2 = size(W,1)/2;  % half-length of spike waveform in samples

% pre-compute the convolution of the waveforms with themselves
[wProj,wNorm] = compWprojW(W);

% load first block of electrode data
ywin = [0,nperchnk];
Y = loadYfun(ywin);

% Run BP on first chunk
ii = ywin(1)+1:ywin(2); % indices to use
Xhat = X0*0; % Initialize Xhat to zeros
[Xhat(ii,:),nrefviols] = runBinaryPursuit(X0(ii,:),Y,W,pspike,wProj,wNorm,minISI);
    
for ichunk = 2:nchnks

    % Report progress (if desired)
    if mod(ichunk,verbose)==0
        fprintf('Estimating spikes (chunk %d of %d)\n', ichunk, nchnks);

    end

    % Determine indices for this chunk of data
    i1 = (ichunk-1)*nperchnk; % index of first bin in this chunk
    i2 = min(ichunk*nperchnk,slen);

    % Load waveforms, if necessary
    wnumneeded = ceil(i2/nperblck);
    if wnum ~= wnumneeded
        wnum = wnumneeded;
        load(sprintf(wfile, wnum),'W');
        [wProj,wNorm] = compWprojW(W); % pre-compute waveform convolutions
    end

    % load electrode data
    ywin = [i1-nw2,i2]; % indices to process
    Y = loadYfun(ywin);

    % perform BP sorting
    [xsrt,nref] = runBinaryPursuit(X0(ywin(1)+1:ywin(2),:),Y,W,pspike,wProj,wNorm,minISI);
    Xhat(i1+1:i2,:) = xsrt(nw2+1:end,:); % insert chunk into Xhat
    nrefviols = nrefviols+nref; % count if any refractory period violations removed
    
end


