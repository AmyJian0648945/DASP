
function [TDOAest] =TDOA_corr_separateSource()
%%% DESCRIPTION:
%	Calculates the Cross-correlation-based time-difference of arrival (TDOA) estimation


computed_rir = load('Computed_RIRs.mat');
numOfMicrophones = length(computed_rir.m_pos);

%% Computing sample delay between the direct path components of the two RIRs.
% This computation is done by indentifying index of first element that is
% =/= from 0 for each RIR and then subtracting between adjacent microphone 
numOfSources = length(computed_rir.s_pos);
index = ones(numOfSources,numOfMicrophones);
delay = ones(numOfSources,numOfMicrophones-1);
for k=1:1:numOfSources
	for j=1:1:numOfMicrophones
		i = 1;
		while (computed_rir.RIR_sources(i,j,k) <= 0) i = i+1; end
		index(k,j) = i;
	end
	
	for j=1:1:numOfMicrophones-1
		delay(k,j) = index(k,j+1)-index(k,j);
	end
end


% TDOAgndTruth = sample delay between the direct path components of the two
% RIRs = ground truth
TDOAgndTruth = delay';

% Generates the Microphone signal
[micSource micNoise] = computeMicSig_separateSources(computed_rir);

%% Estimation of TDOA via cross-correlation
% Step 1: Cross Correlation
r = [];
lag = [];
for i=1:1:numOfSources
	[r_1, lag_1] = xcorr(micSource(:,i, 1), micSource(:,i,2));
	r = [r ;r_1'];
	lag = [lag ;lag_1];
	
end 

r=r'; % Each source is its own column 
lag=lag';

% Step 2: Peak Detection
peakValue = ones(1,numOfSources);
peakLocation = ones(1,numOfSources);
TDOAest = ones(1,numOfSources);
TDOAestError = ones(1,numOfSources);
for i=1:1:numOfSources
	[peakValue(i), peakLocation(i)] = max(r(:,i));
	% Step 3: Locate the delay
	TDOAest(i) = lag(peakLocation(i), i);
	% Chec=k the estimation against the ground truth 
	TDOAestError(i) = TDOAest(i) - TDOAgndTruth(i);
end


% % Plot the Cross-correlation function, with ground truth marked
% figure
% hold on
% plot(r)
% % stem(find(lag == TDOAgndTruth), peakValue) 
% hold off
end
