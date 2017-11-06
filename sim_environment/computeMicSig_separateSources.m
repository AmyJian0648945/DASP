%%% DESCRIPTION:
%	Input: Computed_RIRs.mat (GUI .mat file)
% 	Output: micSource(row=signal, column=the different source, z=different mics)
% 
%	Note: Amount of noise sources in the GUI should be the same as the
%	number of .wav noise sources - if there are multiple GUI noise sources,
%	load the .wav file multiple times


function [micSource, micNoise] = computeMicSig_separateSources(computed_rir)

% Initialisation
numOfMics = size(computed_rir.RIR_sources,2);
numOfSources = size(computed_rir.RIR_sources,3);
numOfNoiseSources = size(computed_rir.v_pos,1); 
source_speech = cell(numOfSources,2);
source_noise = cell(numOfNoiseSources,2);

% User defined variables (!!!)
lenMicSig = 5;		% length of desired microphone signal [sec]

% User defined noise and speech source (!!!)
for i=1:1:numOfSources
	[source_speech{i,1},source_speech{i,2}] = audioread('speech2.wav');

	% [source_speech{2,1},source_speech{2,2}] = audioread('speech2.wav');
	% [source_noise{1,1},source_noise{1,2}] = audioread('White_noise1.wav');
	% [source_noise{2,1},source_noise{2,2}] = audioread('Babble_noise1.wav');
end

samplesToKeep = computed_rir.fs_RIR.*lenMicSig;
for i = 1:numOfSources
	source_speech{i,1} = resample(source_speech{i,1}, computed_rir.fs_RIR, source_speech{i,2});
	source_speech{i,1} = source_speech{i,1}(1:samplesToKeep);
end

for i = 1:numOfNoiseSources
	source_noise{i,1} = resample(source_noise{i,1},computed_rir.fs_RIR,source_noise{i,2});
	source_noise{i,1} = source_noise{i,1}(1:samplesToKeep);
end

micSource = []; micNoise = [];
leng = length(conv(source_speech{1,1}, computed_rir.RIR_sources(:,1,1)));
for i = 1:numOfMics
	if(numOfSources > 0)
		micSource(:,1,i) = conv(source_speech{1,1}, computed_rir.RIR_sources(:,i,1));
	end
	if(numOfSources > 1)
		for j = 2:1:numOfSources 
			micSource(:,j,i) = conv(source_speech{1,1}, computed_rir.RIR_sources(:,i,j));
		end
	end 
	
	if(numOfNoiseSources > 0)
		micNoise(:,1,i) = conv(source_noise{1,1}, computed_rir.RIR_noise(:,i,1));
	end
	if(numOfNoiseSources > 1)
		for k = 2:1:numOfNoiseSources
			micNoise(:,k,i) = conv(source_noise{1,1}, computed_rir.RIR_noise(:,i,k));
		end
	end 
	
end

end