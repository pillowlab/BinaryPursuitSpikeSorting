% Generate some fake data to illustrate the use of binary pursuit spike sorting code.

fprintf('script0_simulateDataForTesting: generating dataset for spike sorting...\n');

% Set params for simulated dataset using those in the 'setSpikeSortParams.m' script
setSpikeSortParams;  % (LOOK HERE FOR DETAILS).

% Set parameters governing the simulated spike trains
nwt = 30;  % # time samples in spike waveforms
sprate = 100; % mean spike rate (unrealistically high to observe lots of simultaneous spks)
nsesig = .1; % marginal stdv of additive noise (arbitrary units)


%% 2.  Generate spike trains % ----------------------
Xsp = double(sparse((rand(sdat.nsamps,sdat.ncell) < sprate/sdat.samprate))); % each column is spk train

% remove spikes too close together
minisi = nwt*1.5; % smallest allowed interspike interval
for j = 1:sdat.ncell
    tsp = find(Xsp(:,j)); % spike times
    isi = [nwt; diff(tsp)]; % interspike intervals
    kk = find(isi<minisi); 
    Xsp(tsp(kk),j) = 0;  % remove these spikes
end


%% 3. Generate some Spike waveforms % ----------------------
tbins = (1:nwt)';  % time bins

% Make some basis waveforms 
someWaves = normpdf(repmat(tbins,1,5),repmat(nwt/5+(1:2:9),nwt,1),repmat(1.5:6,nwt,1));

% Make waveform for each neuron
W = zeros(nwt,sdat.ne,sdat.ncell); % tensor for spike waveforms (nwt x nelectrodes x ncells)
for j = 1:sdat.ncell
    elecinds = max(1,round(j/sdat.ncell*sdat.ne)-2):min(sdat.ne,round(j/sdat.ncell*sdat.ne)+1); % which electrodes the cell talks to
    spwaveform = someWaves*(randn(size(someWaves,2),length(elecinds))+.1); % the waveform
    W(:,elecinds,j) = spwaveform./norm(spwaveform(:)); % put unit-vector waveform in tensor
end

% ----------- MAKE FIG -------------
% Plot all waveforms as [time x electrode] image
subplot(121);
imagesc(reshape(permute(W,[1 3 2]),nwt*sdat.ncell,[])); % image showing all waveforms
set(gca,'ytick',nwt/2:nwt:nwt*sdat.ncell,'yticklabel',1:sdat.ncell);
title('spike waveforms for each neuron');
xlabel('electrode #'); ylabel('cell');
% Plot waveforms for each neuron
for j = 1:sdat.ncell
    subplot(sdat.ncell,2,j*2);
    plot(W(:,:,j));
    set(gca,'xtick', []);
    ylabel(sprintf('cell %d',j));
    if j ==sdat.ncell
        xlabel('time');
    end    
end

%% 4. Simulate recorded electrode data % ----------------------

% Compute noiseless electrode signal by convolving spike trains with waveforms
y0 = compVpredictionSprse(Xsp,W); 

% Make noise filter (for adding realistic noise)
nsefilter_t = exp(-(0:9))';
nsefilter_x = normpdf(-sdat.ne/2:sdat.ne/2,0,1)./normpdf(0);
nsefilter = nsefilter_t*nsefilter_x;
nsefilter = nsefilter./norm(nsefilter(:));

addednoise = conv2(randn(sdat.nsamps,sdat.ne),nsefilter,'same')*nsesig;  % make colored noise
Y = y0+addednoise; % add noise

% ----------- MAKE FIG -------------
% Plot noiseless (red) and noisy (blue) traces for a few electrodes
T = 500;  % time range to display
iipl = 1:T; % time indices to plot
npl = min(10,sdat.ne+1);
for j = 1:npl-1
    subplot(npl,1,j+1);
    plot(iipl, Y(iipl,j), iipl, y0(iipl,j), 'r');
    ylabel(sprintf('electr #%d', j));
    if j < npl-1
        set(gca, 'xticklabel', []);
    end
    box off;
end
xlabel('time (samples)');
% --- Plot spikes ----
subplot(npl,1,1); 
cla; hold on;
for j = 1:sdat.ncell
    tsp = find(Xsp(1:500,j));
    plot([0 T],j*[1 1],'k--');
    if ~isempty(tsp)
        plot(tsp,j,'b.');
    end
end
hold off;
set(gca,'ydir','reverse','ytick',1:sdat.ncell,'ylim', [0 sdat.ncell+1]);
title('spike trains'); ylabel('neuron');

%% 5. Create "initialization" spike train with simultaneous spikes removed % ----------------------

shortisi = nwt/3; % define what counts as "near-simultaneous"
XspTot = sum(Xsp,2); % aggregate spike train
tspTot = find(XspTot); % spike times of aggregate spike triain
isiTot = [shortisi;diff(tspTot)]; % isis 
shortIsis = find((isiTot<shortisi)); % indices of 2nd spike in short isis
shortIsis = setdiff(union(shortIsis-1,shortIsis),0); % all indices

% Create "initial" spike train with partial spikes present
Xsp_init = Xsp;
Xsp_init(tspTot(shortIsis),:) = 0;  % remove these spikes

% Femove 10% of additional spikes
Xsp_init = Xsp_init.*(rand(sdat.nsamps,sdat.ncell)<0.9);

fprintf('----\nTotal number of true spikes: %d\n', full(sum(Xsp(:))));
fprintf('Number of near-synchronous spikes: %d\n', full(sum(XspTot(tspTot(shortIsis)))));
fprintf('Number total spikes removed: %d\n', full(sum(Xsp(:)-Xsp_init(:))));
fprintf('Number spikes in ''Xsp_init'': %d\n', full(sum(Xsp_init(:))));


%% 6. Save out data  % ----------------------
if ~exist('dat','dir');
    mkdir('dat');
end
dirname = sprintf('dat/simdata');
if ~exist(dirname,'dir')
    mkdir(dirname);
end

save([dirname,'/Xsp_true.mat'],'Xsp');
save([dirname,'/W_true.mat'],'W');
save([dirname,'/Y.mat'],'Y');
save([dirname,'/Xsp_init.mat'],'Xsp_init');
