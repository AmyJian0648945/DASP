%%% DESCRIPTION:
%	Input:
%	- Computed_RIRs.mat (GUI .mat file)
% 	- length of the mic signal
%	- flag_output (shows how microphone signal should be outputted)
%		flag_output == 1: Output is [mic,~,~] = Source + Noise;
%		flag_output == 2: Output is [~, micSource, micNoise]
%		flag_output == 3: Output is [~, micSource, micNoise] + 10% power of the power 
%	- flag_input (shows how source / noise should be formatted in the output)
%		flag_input  == 1: Input sources are all name_source1
%		flag_input  == 2: Input source 1 = name_source1, input source 2 = name_source2
%		flag_input  == 3: input noise are all name_noise1
%		flag_input  == 4: Input noise 1 = name_noise1, input noise 2 = name_noise2
% 	Output: mic (variable; conv(noise, noiseRIR) + conv(sig, sigRIR)
% 
%	Note: Amount of noise sources in the GUI should be the same as the
%	number of .wav noise sources - if there are multiple GUI noise sources,
%	load the .wav file multiple times
% 
%	Note: The method of truncation that we used might yield an error


function [mic, micSource, micNoise] = computeMicSig(computed_rir,lenMicSig, flag_output, flag_input, name_source1, name_source2, name_noise1, name_noise2)

% Initialisation
numOfMics = size(computed_rir.RIR_sources,2); %% number of microphones in the scenario
numOfSources = size(computed_rir.s_pos,1); %% number of signal sources in the scenario
numOfNoiseSources = size(computed_rir.v_pos,1); %% number of noise sources in the scenario
source_speech = cell(numOfSources,2); % defining cells for audio sources 
source_noise = cell(numOfNoiseSources,2);% defining cells for noise sources 

% User defined noise and speech source (!!!)In this case every audio sources
% is playing 'speech2.wav' and every noise source is playing 'White_noise1.wav'
% source_speech{i,1} is the signal itself while source speech is
% source_speech{i,2} is the sampling frequency fs

for i=1:1:numOfSources
	[source_speech{i,1},source_speech{i,2}] = audioread('speech1.wav');
end
% [source_speech{1,1},source_speech{1,2}] = audioread('speech1.wav');
% [source_speech{2,1},source_speech{2,2}] = audioread('speech2.wav');

for i=1:1:numOfNoiseSources
	[source_noise{i,1},source_noise{i,2}] = audioread('Babble_noise1.wav');
end
[source_speech{2,1},source_speech{2,2}] = audioread('speech2.wav');
% resampling the signal at fs = fs_RIR instead of the sampling frequency of
% the source speech + truncating the signal (given duration in lenMicSig )
samplesToKeep = computed_rir.fs_RIR.*lenMicSig;
for i = 1:1:numOfSources
	source_speech{i,1} = resample(source_speech{i,1},computed_rir.fs_RIR,source_speech{i,2}); %% resampling
	source_speech{i,1} = source_speech{i,1}(1:samplesToKeep); %% truncation 
end

for i = 1:1:numOfNoiseSources
	source_noise{i,1} = resample(source_noise{i,1},computed_rir.fs_RIR,source_noise{i,2});%% resampling
	source_noise{i,1} = source_noise{i,1}(1:samplesToKeep); %% truncation 
end


leng = length(filter(computed_rir.RIR_sources(:,1,1),1,source_speech{1,1}));
mic = zeros(leng,numOfMics); % preallocation of mic siganls-> each column is a mic signal and there's
% numOfMics columns 
micSource = []; micNoise = [];


if(flag_output == 1 || flag_output == 3)
	for i = 1:1:numOfMics
		speech = zeros(leng,1);
		if(numOfSources > 0)
			speech = filter(computed_rir.RIR_sources(:,i,1),1,source_speech{1,1}); %% audio signal is filtered by the Room IR
		end
		if(numOfSources > 1)
			for j = 2:1:numOfSources 
				speech = speech + filter(computed_rir.RIR_sources(:,i,j),1,source_speech{j,1}); % each mic received a weighted sum of audio signal. Weights are the RIR 
			end
		end 

		noise = zeros(leng,1);
		if(numOfNoiseSources > 0)
			noise = filter(computed_rir.RIR_noise(:,i,1),1,source_noise{1,1}); %% noise signal is filtered by the Room IR
		end
		if(numOfNoiseSources > 1)
			for k = 2:1:numOfNoiseSources
				noise = noise + filter(computed_rir.RIR_noise(:,i,k),1,source_noise{k,1}); % each mic received a weighted sum of noise signal. Weights are the RIR
			end
		end 

		if(flag_output==1) 
			mic(:,i) = speech + noise; % total signal received by the mic is AUDIO + NOISE
		end
	end
	
	if(flag_output == 3) %%% For Session 3 %%%
		
		% Identify where speech signals are active
		VAD = abs(speech(:,1))>std(speech(:,1))*1e-3;

		% Find the power of the signal of mic1 (where its active)
		speechPower = var(speech(VAD==1,1));
		noiseVar = 0.1*speechPower;
		addWhiteGausNoise = random('Normal', 0, sqrt(noiseVar), size(speech,1), size(speech,2));

		% Add the noise to 'speech'
		mic = speech + addWhiteGausNoise + noise; 

		% Compute SNR
		noisePower = var(addWhiteGausNoise(:,1) + noise(:,1));
		SNR = 10*log10(speechPower ./ noisePower)

		save('SNR_in','SNR');
		save('mic','mic','speech','addWhiteGausNoise','SNR');
		micSource = speech; 
		micNoise = addWhiteGausNoise + noise; 
	end
		
		
elseif(flag_output == 2)
	
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

	

else 
	error('Somethings wrong with the flag_output')
end
	
	
end