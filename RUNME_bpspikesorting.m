% Example script to generate a simulated dataset and then illustrate the steps in our
% spike sorting algorithm.

% 0. Generate a simulated dataset
script0_simulateDataForTesting;

%% 1. Estimate spike waveforms for each neuron using spike times provided
step1_estimWaveforms

%% 2. Estimate the temporal and spatial noise covariance and whiten the raw data
step2_WhitenData;

%% 3. Re-estimate spike waveforms using whitened data
step3_reestimWaveforms;

%% 4. Run binary pursuit: identify the spike times given waveforms and whitened data
step4_BPspikesort;

%% 5. Compare simulated and estimated spike trains (ONLY RELEVANT FOR SIMULATED DATA)
step5_analyzePerformance_simData;


% ===============================================================================
% NOTE: To run this on your own data, you'll need to:
%
% (1) specify the number of neurons and identify
% enough spike times from each of them to get a reasonable first-pass estimate of the
% spike waveforms (e.g., using clustering).
%
% (2) Specify a few dataset and processing dependent parameters (e.g., number of
% electrodes, number of neurons, sample rate, number of time bins in spike waveforms, 
% path to the raw data, a function for loading time windowed chunks of data) in the
% general script 'setSpikeSortParams.m' 
%
% I recommend opening up each of the step_x scripts above to get a sense of what's
% happening in each step.  Please feel free to email me with comments, questions,
% suggestions and bug reports (pillow@princeton.edu).
% ===============================================================================
