%  step2_whitenData.m
%
% Compute whitening filters with current residuals and apply to raw electrode data

% Set path and loads relevant data structures: 'sdat','dirlist', 'filelist'
setSpikeSortParams;  

% --- Load initial estimate of spike times (sparse nsamps x ncell array) ----------
Xsp = struct2array(load(filelist.initspikes));  % loads variable 'Xsp_init'

% ---  Set params governing block size ------
nsamps = sdat.nsamps; % Number of total samples
nsampsPerW = sdat.nsampsPerW;  % number of samples to use for each estimate
nWblocks = ceil(nsamps/nsampsPerW);  % number of blocks for estimating waveforms

% --- Set params governing the whitening ------
nxc_t = 16;  % # time bins for temporal whitening filter
nxc_x = 5;  % # time bins to use while spatial whitening (not too big, <= 9)

% --- Determine number of "chunks" to divide whitened data into (for BP sort) -----
nsPerChunk = sdat.nsampsPerBPchunk;  % number of samples per chunk
nchunksPerBlock = nsampsPerW/sdat.nsampsPerBPchunk; % number of chunks per "waveform block"


% Loop over blocks
for blocknum = 1:nWblocks
    fprintf('Step 2: whitening residuals (block %d/%d)\n',blocknum,nWblocks);
    
    % --- Set time window and load relevant data ----
    twin = [(blocknum-1)*nsampsPerW, min(blocknum*nsampsPerW,nsamps)]; % time window
    Ydat = loadrawY(twin); % load the relevant block of electrode data
    load(sprintf(filelist.Wraw, blocknum),'W'); % load spike waveform
    
    %---  Compute whitening filters & whiten electrode data ------------------
    [Ywht,tfilts,xfilts] = ...
        compWhitening(Xsp(twin(1)+1:twin(2),:),Ydat,W,nxc_t,nxc_x);

    %  Save out whitened data in chunks
    for ichunknum = 1:nchunksPerBlock;
        chunknum = (blocknum-1)*nchunksPerBlock+ichunknum;
        savename = sprintf(filelist.Ywht, chunknum);
        Y = Ywht(((ichunknum-1)*nsPerChunk+1):(ichunknum*nsPerChunk),:);
        save(savename, 'Y');
    end

end
