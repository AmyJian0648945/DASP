clear all; 
%%% DESCRIPTION:
%	Calculates the Cross-correlation-based time-difference of arrival (TDOA) estimation


computed_rir = load('Computed_RIRs.mat');

% Finds first point that is not zero
i = 1;
while (computed_rir.RIR_sources(i,1) <= 0) i = i+1; end
first_index = i;
i = 1;
while (computed_rir.RIR_sources(i,2) <= 0) i = i+1; end

% TDOAgndTruth = sample delay between the direct path components of the two
% RIRs = ground truth
TDOAgndTruth = abs(i-first_index);

% Generates the Microphone signal
mic = computeMicSig(computed_rir);

%% Estimation of TDOA via cross-correlation
% Step 1: Cross Correlation
[r, lag] = xcorr(mic(:,1), mic(:,2));

% Step 2: Peak Detection
[~, peakLocation] = max(r);

% Step 3: Locate the delay
TDOAest = abs(lag(peakLocation))


% Check the estimation against the ground truth 
TDOAestError = TDOAest - TDOAgndTruth

