%  step4_BPspikesort.m
%  --------------------
%  Estimate spike times on whitened data using Binary Pursuit

% Set path and loads relevant data structures: 'sdat', 'dirlist', 'filelist'
setSpikeSortParams;  
fprintf('Step 4: running Binary Pursuit to estimate spike times\n');

% ---- Load initial estimate of spike times (sparse nsamps x ncells array) ------------
X0 = struct2array(load(filelist.initspikes));  % loads variable 'Xsp_init'

% ---  Set some params governing block size for estimating waveforms ------
nsamps = sdat.nsamps; % Number of total samples
nsampsPerW = sdat.nsampsPerW;  % number of samples to use for each estimate
nWblocks = ceil(nsamps/nsampsPerW);  % number of blocks for estimating waveforms
blksize = sdat.nsampsPerBPchunk; % number of samples to process at once for BP. (default 10K)

% Set prior probability of a spike for each neuron (using base rate in initial sort)
pspike0 = mean(X0); % prior probability of a spike for each neuron 

% -- Do sorting ----
[Xhat,nrefviols] = estimSps_BinaryPursuit(loadwhitenedY,filelist.Wwht,...
    nsampsPerW,X0,pspike0,blksize,sdat.minISIsamps);

% -- Report refractory period violations (if necessary) ----
if sum(nrefviols) == 0
    fprintf('No refractory period violations (%.1fms)\n',sdat.minISIms);
else
    fprintf('Number spikes pruned due to refractory violations (%.1fms)\n',sdat.minISIms);
    for j = 1:sdat.ncells
        if nrefviols(j)>0
            fprintf('cell %d: %d\n', j, nrefviols(j));
        end
    end
end

% -- Save out ------
save(filelist.Xhat,'Xhat','nrefviols');

%  ============= NOTES ================
%
% Note 1: If getting too many spikes out (i.e., unrealistically many spikes), try dividing
% pspike0 10, 10^2 10^3, etc.

% Note 2: for an assessment of the reliability of each neuron's spikes (i.e., sorting
% accuracy) , try dividing pspike0 by 2 or multiplying it by 2, and see how many more or
% fewer spikes you get. For a highly reliable sort (i.e., where the posterior is strongly
% determined by the likelihood), you should get nearly the same answer, regardless of the
% prior. 

