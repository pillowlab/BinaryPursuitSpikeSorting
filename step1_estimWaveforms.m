%  step1_estimWaveforms.m
%  -------------------------
%  Estimate spike waveforms from initial subset of spikes found (before whitening)

% Set path and loads relevant data structures: 'sdat', 'dirlist', 'filelist'
setSpikeSortParams;  

% Is simulation?  (for comparisons to ground truth)
isSIM = 1; % Set to true only when running demo code with simulation data 

% ---- Load initial estimate of spike times (sparse nsamps x ncells array) -----
Xsp = struct2array(load(filelist.initspikes));  % loads variable 'Xsp_init'

% ---  Set some params governing block size for estimating waveforms ------
nsamps = sdat.nsamps; % Number of total samples
nsampsPerW = sdat.nsampsPerW;  % number of samples to use for each estimate
nWblocks = ceil(nsamps/nsampsPerW);  % number of blocks for estimating waveforms


%% Estimate waveform independently on each chunk of data (here 30s worth / chunk)
for blocknum = 1:nWblocks
    fprintf('Step 1: estimating pre-whitened waveforms (block %d/%d)\n', blocknum,nWblocks);
    
    % --- Set time window and load electrode data ---
    twin = [(blocknum-1)*nsampsPerW, min(blocknum*nsampsPerW,nsamps)]; % time window
    Ydat = loadrawY(twin); % load the relevant block of electrode data
    
    % --- Estimate spike waveform ----
    W = estimWaveforms(Xsp(twin(1)+1:twin(2),:),Ydat,sdat.nw); 
    
    % --- Save out waveforms ---
    savename = sprintf(filelist.Wraw, blocknum);
    save(savename, 'W', 'twin');

end

%% Make plot of estimated waveforms

% NOTE: if support of waveform extends outside plotted window, consider shifting spike
% times or increase sdat.nsampsPerW
for j = (1:sdat.ncells)
    subplot(sdat.ncells,1,j);
    plot(1:sdat.nw, W(:,:,j),'b');
    ylabel(sprintf('cell %d',j));
end
xlabel('time (bins)');

% === Compare esimated and true waveform (SIMULATION ONLY) === 
if isSIM 
    ww = load('dat/simdata/W_true.mat');  % load true waveforms
    for j = (1:sdat.ncells)
        subplot(sdat.ncells,1,j); hold on;
        plot(1:size(ww.W,1),ww.W(:,:,j),'r--'); hold off;
    end
end