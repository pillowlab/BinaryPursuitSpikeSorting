% PATH SETTING
addpath sort-tools/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMS GOVERNING THE RECORDING DATA (sdat = DAT STRUCT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHANGE THESE FOR YOUR DATASET
sdat.ne = 6; % number of electrodes
sdat.ncell = 4; % number of neurons
sdat.nsecs = 120;  % number of total seconds in the recording
sdat.samprate = 20e3; % sample rate
sdat.nsamps = sdat.nsecs*sdat.samprate;  % number of total samples


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET PARAMS GOVERNING PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CAN LEAVE THESE FIXED, BUT CHANGE TO OPTIMIZE PERFORMANCE
sdat.nsampsPerFile= 20000; % num samples per saved (processed) electrode-data file
sdat.nsecsPerW = 60; % number of seconds' data to use for each estimate of spike waveform
sdat.nsampsPerW = sdat.nsecsPerW*sdat.samprate; % num samples per waveform estimate 
sdat.nw = 30; % number time bins in spike waveform (MUST BE EVEN)
sdat.nsampsPerBPchunk = 10000; % number samples in a chunk for binary pursuit sorting
sdat.minISIms = 0.5; % minimum ISI, in milliseconds
sdat.minISIsamps = round(sdat.minISIms*sdat.samprate/1000);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORKING DIRECTORIES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directory where raw data to be loaded from
dirlist.rawdat = sprintf('dat/simdata/');  % CHANGE THIS FOR YOUR DATASET

% Directories where intermediate processing data to be stored
dirlist.procdat = 'dat/procdat/'; % directory for processed data (to be created)
dirlist.W = [dirlist.procdat 'Wraw/'];  % raw waveform estimates (pre-whitening)
dirlist.Wwht = [dirlist.procdat 'Wwht/'];  % waveform estimates after whitening
dirlist.Ywht = [dirlist.procdat 'Ywht/'];  % sparsified waveform estimates
dirlist.tspEstim = [dirlist.procdat 'tspEstim/'];  % spike train estimates

% --------  Check that all dirs have been created ------------
dirfields = fieldnames(dirlist);
for jj = 1:length(dirfields)
    dirname = dirlist.(dirfields{jj}); % dynamic field names (may break older matlab versions)
    if ~isdir(dirname);
        fprintf('SETSPIKESORTPARAMS: making directory ''%s''\n', dirname);
        mkdir(dirname);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WORKING FILENAMES  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NAMES FOR RAW DAT FILES (CHANGE THESE)
filelist.initspikes = [dirlist.rawdat, 'Xsp_init.mat']; % initial spike train estimate (sparse nsamps x ncells array)
filelist.Ydat = [dirlist.rawdat, 'Y.mat']; % initial spike train estimate (sparse nsamps x ncells array)

% NAMES FOR PROCESSED DATA FILES (can leave)
filelist.Ywht = [dirlist.Ywht, 'Y_chunk%d.mat']; % initial spike train estimate (sparse nsamps x ncells array)
filelist.Wraw = [dirlist.W, 'Wraw_%d.mat']; % initial (pre-whitening) estimates of Waveform 
filelist.Wwht = [dirlist.W, 'Wwht_%d.mat']; % estimates of Waveform w/ whitened data
filelist.Xhat = [dirlist.tspEstim, 'Xhat.mat'];

% Function for loading raw electrode data (WRITE THIS FOR YOUR OWN DATA!)
loadrawY = @(twin)(loadRawElecDatWin(twin,filelist.Ydat));

% Function for loading whitened electrode data (no need to rewrite)
loadwhitenedY = @(twin)(loadWhiteElecDatWin(twin,filelist.Ywht,sdat.nsampsPerBPchunk));

