% step5_analyzePerformance_simData.m
% --------------------
% Script to examine performance on simulated data
%
% Note that you can vary performance on simulated data by changing the 'nsesig' 
% parameter in the script 'script0_simulateDataForTesting.m', which controls SNR

% Set path and loads relevant data structures: 'sdat', 'dirlist', 'filelist'
setSpikeSortParams;  

% ---- Load initial estimate of spike times (sparse nsamps x ncell array) ------------
Xtrue = struct2array(load('dat/simdata/Xsp_true.mat')); % true spikes
Xhat = struct2array(load(filelist.Xhat,'Xhat')); % estimated spikes

% Set tolerance in spike time for counting a hit or misses
tol = 3;  % in bins

% Compute Hits for each time shift
slen = sdat.nsamps;
Xhit = zeros(2*tol+1,sdat.ncell);
XFA = zeros(1,sdat.ncell);
for j = 1:(2*tol+1)
    ii1 = max(1,j-tol):min(slen,slen-tol+j-1);
    ii2 = max(1,tol-j+2):min(slen,slen-j+tol+1);
    Xhit(j,:) = sum(Xtrue(ii1,:).*Xhat(ii2,:));
end

%% Report Quantitative Performace 

% Make fig showing p(Hit) vs. time shift
nHits = sum(Xhit);
FracHits = Xhit./repmat(sum(Xtrue),2*tol+1,1);
clf;
plot(-tol:tol, FracHits);
title('P(Hit) vs. time shift');
xlabel('time shift (bins)'); ylabel('P(hit)');

% Report Misses
nMiss = sum(Xtrue)-nHits;
fracMiss = nMiss./sum(Xtrue);
fprintf('\n--- Misses ----\n');
for j = 1:sdat.ncell
    fprintf('cell %d: %d (frac = %.3f)\n', j,nMiss(j), fracMiss(j));
end

% Report FAs
nFA = sum(Xhat)-sum(Xhit);
fracFA = nFA./sum(Xhat);
fprintf('\n--- False Alarms ----\n');
for j = 1:sdat.ncell
    fprintf('cell %d: %d (frac = %.3f)\n', j,nFA(j), fracFA(j));
end
