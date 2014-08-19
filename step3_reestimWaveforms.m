%  step3_reestimWaveforms.m
%  -------------------------
%  Estimate spike waveforms from initial subset of spikes  (after whitening)

% Set path and loads relevant data structures: 'sdat', 'dirlist', 'filelist'
setSpikeSortParams;  

% ---- Load initial estimate of spike times (sparse nsamps x ncells array) ------------
Xsp = struct2array(load(filelist.initspikes));  % loads variable 'Xsp_init'

% ---  Set some params governing block size for estimating waveforms ------
nsamps = sdat.nsamps; % Number of total samples
nsampsPerW = sdat.nsampsPerW;  % number of samples to use for each estimate
nWblocks = ceil(nsamps/nsampsPerW);  % number of blocks for estimating waveforms


%% Estimate waveform independently on each chunk of data (here 30s worth / chunk)
for blocknum = 1:nWblocks
    fprintf('Step 3: estimating whitened waveforms (block %d/%d)\n', blocknum,nWblocks);
    
    % --- Set time window and load relevant data ----
    twin = [(blocknum-1)*nsampsPerW, min(blocknum*nsampsPerW,nsamps)]; % time window
    Y = loadwhitenedY(twin);
    
    % --- Estimate spike waveform ----
    [W,wsigs] = estimWaveforms(Xsp(twin(1)+1:twin(2),:),Y,sdat.nw); 
    
    % --- Prune waveforms --- 
    % Estimate the electrode support for each waveform 
    % (which electrodes carry non-zero signal)
    % Not implemented for now.
    
    % --- Save out waveforms ---
    savename = sprintf(filelist.Wwht, blocknum);
    save(savename, 'W', 'twin');

end

%% Make plot of estimated waveforms

% NOTE: if support of waveform extends outside plotted window, consider shifting spike
% times or increase sdat.nsampsPerW
for j = (1:sdat.ncells)
    subplot(sdat.ncells,1,j);
    plot(1:sdat.nw, W(:,:,j),'b'); axis tight;
    ylabel(sprintf('cell %d',j)); 
end
xlabel('time (bins)');
