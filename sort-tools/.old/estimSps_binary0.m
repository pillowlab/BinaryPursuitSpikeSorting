function xxEst = estimSps_binary(XX0,yyfile,WWname,pspike,sortwin,nsPerCh);
% xxEst = estimSps_binary(XX0,yyfile,WWname,pspike,sortwin,nsecsPerChunk);
%
%  Estimates spike train over a large chunk of data by calling:
%   - estimSps_binary_chunk  (binary sorting)
%  on smaller chunks of data
%
%  Inputs:  XX0 = initial guess at spike trains
%           yyfile = name of file with electrode data
%           WWname = name of file containing waveforms
%           pspike = prior probability of a spike in a bin
%           sortwin = time window over which to estimate spike train (in bins)
%           nsecsPerChunk = chunk size (in secs) per optimization step
%
% Outputs:  
%   xxEst = estimated binary spike train

global SampRate TWIN_EXPT MINSPERWW;
 
% --- Load initial waveform to get information about ---
wfilename = sprintf(WWname, 1);
if ~exist(wfilename)
    error(sprintf('Waveform filename does not exist: %s', wfilename));
else
    load(wfilename);
end
[nw,ne,nc] = size(WW);  % # time samples, # electrodes, # cells

% Set some up params
tbuff = ceil(nw*2);  % # of bins buffer: discard for each estimated chunk
xlen = diff(sortwin);
xxEst = sparse(xlen,nc);
chunklen = nsPerCh*SampRate;
nchunks = ceil(xlen/chunklen);

for jchunk = 1:nchunks
    fprintf('--- estimSps_binary: chunk %d (%ds) ---\n',jchunk,nsPerCh);
    
    % --- Define time window -----
    xwin = [0 chunklen] + (jchunk-1)*chunklen; % Window to save (output)
    twin = xwin + tbuff*[-1 1]; % Window for electrode dat (input)
    twin = [max(twin(1),TWIN_EXPT(1)), min(twin(2),TWIN_EXPT(2))];
    
    % ---- spike waveforms ----
    wwblocknum = ceil((xwin(1)+1)/(SampRate*60*MINSPERWW));
    load(sprintf(WWname, wwblocknum));
    
    % ---- Load electrode data ----
    yy = loadElecDatWin(twin,yyfile,ne); % load electr data
    xx = XX0(twin(1)+1:twin(2),:); % spiking data for this chunk
    
    % ----  Do optimization ---------------------
    tic; [xsrt,yresid,dLi]=estimSps_binary_chunk(xx,yy,WW,pspike); toc;
    
    %-----  Extract relevant chunk -------
    iikp = xwin(1)-twin(1)+[0 chunklen];
    xxEst(xwin(1)+1:xwin(2),:) = xsrt(iikp(1)+1:iikp(2),:);    
    
end
