%%% DESCRIPTION:
%	Input: Computed_RIRs.mat (GUI .mat file)
% 	Output: micSource(row=signal, column=the different source, z=different mics)
% 
%	Note: Amount of noise sources in the GUI should be the same as the
%	number of .wav noise sources - if there are multiple GUI noise sources,
%	load the .wav file multiple times
%   The goal of this function is to compute the audio source i seen by a
%   given microphone z through the RIR


function [micSource, micNoise] = computeMicSig_separateSources(computed_rir)

% Initialisation
numOfMics = size(computed_rir.RIR_sources,2); %% number of microphones in the scenario
numOfSources = size(computed_rir.s_pos,1); %% number of signal sources in the scenario
numOfNoiseSources = size(computed_rir.v_pos,1); %% number of noise sources in the scenario
source_speech = cell(numOfSources,2); % defining cells for audio sources 
source_noise = cell(numOfNoiseSources,2);% defining cells for noise sources 

% User defined variables (!!!)
lenMicSig = 5;		% length of desired microphone signal [sec]

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

micSource = []; micNoise = [];
leng = length(conv(source_speech{1,1}, computed_rir.RIR_sources(:,1,1)));
for i = 1:numOfMics
	if(numOfSources > 0)
		micSource(:,1,i) = filter(computed_rir.RIR_sources(:,i,1),1,source_speech{1,1});
	end
	if(numOfSources > 1)
		for j = 2:1:numOfSources 
			micSource(:,j,i) = filter(computed_rir.RIR_sources(:,i,j),1,source_speech{j,1});
		end
	end 
	
	if(numOfNoiseSources > 0)
		micNoise(:,1,i) = filter(computed_rir.RIR_noise(:,i,1),1,source_noise{1,1});
	end
	if(numOfNoiseSources > 1)
		for k = 2:1:numOfNoiseSources
			micNoise(:,k,i) = filter(computed_rir.RIR_noise(:,i,k),1,source_noise{k,1});
		end
	end 
	
end

end