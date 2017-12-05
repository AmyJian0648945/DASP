%%% DESCRIPTION:
%	Input: Computed_RIRs.mat (GUI .mat file)
% 	Output: mic (variable; conv(noise, noiseRIR) + conv(sig, sigRIR)
% 
%	Note: Amount of noise sources in the GUI should be the same as the
%	number of .wav noise sources - if there are multiple GUI noise sources,
%	load the .wav file multiple times
% 
%	Note: The method of truncation that we used might yield an error


function [mic] = computeMicSig(computed_rir,lengMicSig)

% Initialisation
numOfMics = size(computed_rir.RIR_sources,2); %% number of microphones in the scenario
numOfSources = size(computed_rir.s_pos,1); %% number of signal sources in the scenario
numOfNoiseSources = size(computed_rir.v_pos,1); %% number of noise sources in the scenario
source_speech = cell(numOfSources,2); % defining cells for audio sources 
source_noise = cell(numOfNoiseSources,2);% defining cells for noise sources 

% User defined variables (!!!)
lenMicSig = lengMicSig;		% length of desired microphone signal [sec]

% User defined noise and speech source (!!!)In this cas every audio sources
% is playing 'speech2.wav' and every noise source is playing 'White_noise1.wav'
% source_speech{i,1} is the signal itself while source speech is
% source_speech{i,2} is the sampling frequency fs
for i=1:1:numOfSources
	[source_speech{i,1},source_speech{i,2}] = audioread('White_noise1.wav');
end
for i=1:1:numOfNoiseSources
	[source_noise{i,1},source_noise{i,2}] = audioread('White_noise1.wav');
end
[source_speech{2,1},source_speech{2,2}] = audioread('White_noise1.wav');
% resampling the signal at fs = fs_RIR instead of the sampling frequency of
% the source speech + truncating the signal (given duration in lenMicSig )
samplesToKeep = computed_rir.fs_RIR.*lenMicSig;
for i = 1:1:numOfSources
	source_speech{i,1} = resample(source_speech{i,1},computed_rir.fs_RIR,source_speech{i,2}); %% resampling
	source_speech{i,1} = source_speech{i,1}(1:samplesToKeep); %% truncation 
end

for i = 1:1:numOfNoiseSources
	source_noise{i,1} = resample(source_noise{i,1},computed_rir.fs_RIR,source_noise{i,2});%% resampling
	source_noise{i,1} = source_noise{i,1}(1:samplesToKeep);%% truncation 
end


leng = length(filter(computed_rir.RIR_sources(:,1,1),1,source_speech{1,1}));
mic = zeros(leng,numOfMics); % preallocation of mic siganls-> each column is a mic signal and there's
% numOfMics columns 
for i = 1:1:numOfMics
	tempSource = zeros(leng,1);
	if(numOfSources > 0)
		tempSource = filter(computed_rir.RIR_sources(:,i,1),1,source_speech{1,1}); %% audio signal is filtered by the Room IR
	end
	if(numOfSources > 1)
		for j = 2:1:numOfSources 
			tempSource = tempSource + filter(computed_rir.RIR_sources(:,i,j),1,source_speech{j,1}); % each mic received a weighted sum of audio signal. Weights are the RIR 
		end
    end 
    
	tempNoise = zeros(leng,1);
	if(numOfNoiseSources > 0)
		tempNoise = filter(computed_rir.RIR_noise(:,i,1),1,source_noise{1,1}); %% noise signal is filtered by the Room IR
	end
	if(numOfNoiseSources > 1)
		for k = 2:1:numOfNoiseSources
			tempNoise = tempNoise + filter(computed_rir.RIR_noise(:,i,k),1,source_noise{k,1}); % each mic received a weighted sum of noise signal. Weights are the RIR
		end
	end 
	
	mic(:,i) = tempSource + tempNoise; % total signal received by the mic is AUDIO + NOISE
	
end

end