%%% DESCRIPTION:
%	Input:
%	- Computed_RIRs.mat (GUI .mat file)
% 	- length of the mic signal
%	- flag_output (shows how microphone signal should be outputted)
%		flag_output == 1: Output is [mic,~,~] = Source + Noise;
%		flag_output == 2: Output is [~, micSource, micNoise]
%		flag_output == 3: Output is [~, micSource, micNoise] + 10% power of the power 
%	- flag_input (shows how source / noise should be formatted in the output)
%		flag_input  == 1: Sources, Noises are all the same 
%		flag_input  == 2: Noises are all the same, sources are all different
%		flag_input  == 3: Sources are all the same, noises are all different
%		flag_input  == 4: Sources, Noises are all different
% 	Output: mic (variable; conv(noise, noiseRIR) + conv(sig, sigRIR)
% 
%	Note: Amount of noise sources in the GUI should be the same as the
%	number of .wav noise sources - if there are multiple GUI noise sources,
%	load the .wav file multiple times
% 
%	Note: The method of truncation that we used might yield an error


function [mic, micSource, micNoise] = computeMicSig(computed_rir,lenMicSig, flag_output, flag_input, filename_source, filename_noise)

% Testing the sampling frequency --> should be 44.1kHz for session 3 
% if (computed_rir.fs_RIR ~= 44.1e3)
%     error('The sampling frequency is not 44.1 kHz');
% end
% Initialisation
numOfMics = size(computed_rir.RIR_sources,2); %% number of microphones in the scenario
numOfSources = size(computed_rir.s_pos,1); %% number of signal sources in the scenario
numOfNoiseSources = size(computed_rir.v_pos,1); %% number of noise sources in the scenario
source_speech = cell(numOfSources,2); % defining cells for audio sources 
source_noise = cell(numOfNoiseSources,2);% defining cells for noise sources 


% source_speech{i,1} is the signal itself while source speech is
% source_speech{i,2} is the sampling frequency fs

% Formating the sources based on 'flag_input'
if(flag_input == 1) %% All sources are the same
	%%% Sources %%%
	for i=1:1:numOfSources
		[source_speech{i,1},source_speech{i,2}] = audioread(filename_source{1});
	end
	
	%%% Noises %%%
	for i=1:1:numOfNoiseSources
		[source_noise{i,1},source_noise{i,2}] = audioread(filename_noise{1});
	end
elseif(flag_input == 2) %% Sources are different, noises are the same
	
	%%% Sources %%%
	if(numOfSources ~= size(filename_source,2)) 
		error('Number of sources in the GUI must be the same as the number of .wav files given by the user!')
	end
	for i=1:1:numOfSources
		[source_speech{i,1},source_speech{i,2}] = audioread(filename_source{i});
	end
	
	%%% Noises %%%
	for i=1:1:numOfNoiseSources
		[source_noise{i,1},source_noise{i,2}] = audioread(filename_noise{1});
	end
	
elseif(flag_input == 3) %% Sources are the same, Noises are different
	
	%%% Sources %%%
	for i=1:1:numOfSources
		[source_speech{i,1},source_speech{i,2}] = audioread(filename_source{1});
	end
	
	%%% Noises %%%
	if(numOfNoiseSources ~= size(filename_noise,2)) 
		error('Number of sources in the GUI must be the same as the number of .wav files given by the user!')
	end
	for i=1:1:numOfNoiseSources
		[source_noise{i,1},source_noise{i,2}] = audioread(filename_noise{i});
	end
	
elseif(flag_input == 4) %% All sources are different!
	
	%%% Sources %%%
	if(numOfSources ~= size(filename_source,2)) 
		error('Number of sources in the GUI must be the same as the number of .wav files given by the user!')
	end
	for i=1:1:numOfSources
		[source_speech{i,1},source_speech{i,2}] = audioread(filename_source{i});
	end
	
	%%% Noises %%%
	if(numOfNoiseSources ~= size(filename_noise,2)) 
		error('Number of sources in the GUI must be the same as the number of .wav files given by the user!')
	end
	for i=1:1:numOfNoiseSources
		[source_noise{i,1},source_noise{i,2}] = audioread(filename_noise{i});
	end

else
	error('Please enter a flag_input value [1,4]')
end 









% [source_speech{2,1},source_speech{2,2}] = audioread('speech2.wav');
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
mic = zeros(leng,numOfMics);
noiseTot= zeros(leng,numOfMics);
speechTot = zeros(leng,numOfMics);% preallocation of mic siganls-> each column is a mic signal and there's
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
        if(flag_output==3) 
			speechTot(:,i) = speech;
            noiseTot(:,i) = noise;% total signal received by the mic is AUDIO + NOISE
        end
	end
	
	if(flag_output == 3) %%% For Session 3 %%%
		
		% Identify where speech signals are active
        % This compute the standard deviation and take 0.1% of it 
        % 
		VAD = abs(speechTot(:,1))>std(speechTot(:,1))*1e-3;

		% Find the power of the signal of mic1 (where its active)
		speechPower = var(speechTot(VAD==1,1));
		noiseVar = 0.1*speechPower;
		addWhiteGausNoise = random('Normal', 0, sqrt(noiseVar), size(mic,1), size(mic,2));

		% Add the noise to 'speech'
		mic = speechTot+ addWhiteGausNoise + noiseTot; 

		% Compute SNR
		noisePower = var(addWhiteGausNoise(:,1) + noiseTot(:,1));
		SNR = 10*log10(speechPower ./ noisePower);

		save('SNR_in','SNR');
		save('mic','mic','speech','addWhiteGausNoise','SNR');
		micSource = speechTot; 
		micNoise = addWhiteGausNoise + noiseTot; 
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