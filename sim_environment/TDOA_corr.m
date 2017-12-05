function [TDOAest] =TDOA_corr()
%%% DESCRIPTION:
%	Calculates the Cross-correlation-based time-difference of arrival (TDOA) estimation


computed_rir = load('Computed_RIRs.mat');
numOfMicrophones = length(computed_rir.m_pos);


%% Computing sample delay between the direct path components of the two RIRs.
% This computation is done by indentifying index of first element that is
% =/= from 0 for each RIR and then subtracting between adjacent microphone 
index = [];
for j=1:1:numOfMicrophones
	i = 1;
	while (computed_rir.RIR_sources(i,j) <= 0) i = i+1; end
	index = [index i];
end
delay = [];
for j=1:1:numOfMicrophones-1
	delay = [delay index(j+1)-index(j)];
end

% TDOAgndTruth = sample delay between the direct path components of the two
% RIRs = ground truth. Here we have chosen the delay between mic 1 and 2 as
% TDOA ground thruth
TDOAgndTruth = delay(1);

% Generates the Microphone signal
mic = computeMicSig(computed_rir,10);

%% Estimation of TDOA via cross-correlation
% Step 1: Cross Correlation
[r, lag] = xcorr(mic(:,1), mic(:,2));

% Step 2: Peak Detection
[peakValue, peakLocation] = max(r);

% Step 3: Locate the delay
TDOAest = lag(peakLocation);

% Check the estimation against the ground truth 
TDOAestError = abs(TDOAest) - abs(TDOAgndTruth)

% Plot the Cross-correlation function, with ground truth marked
figure
hold on
plot(r)
% stem(find(lag == TDOAgndTruth), peakValue) 
hold off
end
