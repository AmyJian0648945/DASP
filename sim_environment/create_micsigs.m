clear all; 


% Purpose of this file: to read input files, then truncate the file by the user specified lenghth, 
% and resample the signal according to the sampling frequency of the RIR

computed_rir = load('Computed_RIRs.mat'); % always call this before you start, because you need the sampling frequency of the RIR

% User defined variables
lenMicSig = 5; % length of desired microphone signal in seconds
fs_source = [0 0 0 0];
numOfMicsGUI = size(computed_rir.RIR_sources,2);
numOfSourcesGUI = size(computed_rir.RIR_sources,3);
numOfNoiseSourcesGUI = size(computed_rir.v_pos,1);
numInputSource = 1; 
numInputNoise = 0;

source_speech = cell(numInputSource);
source_noise = cell(numInputNoise);

[source_speech{1,1},source_speech{1,2}] = audioread('speech2.wav');
% [source_speech{2,1},source_speech{2,2}] = audioread('speech2.wav');
[source_noise{1,1},source_noise{1,2}] = audioread('Babble_noise1.wav');
% [source_noise{2,1},source_noise{2,2}] = audioread('White_noise1.wav');


samplesToKeep = computed_rir.fs_RIR.*lenMicSig;
for i = 1:numInputSource
	source_speech{i,1} = resample(source_speech{i,1},computed_rir.fs_RIR,source_speech{i,2});
	source_speech{i,1} = source_speech{i,1}(1:samplesToKeep);
end

for i = 1:numInputNoise
	source_noise{i,1} = resample(source_noise{i,1},computed_rir.fs_RIR,source_noise{i,2});
	source_noise{i,1} = source_noise{i,1}(1:samplesToKeep);
end


leng = length(conv(source_speech{1,1}, computed_rir.RIR_sources(:,1,1)));
for i = 1:numOfMicsGUI
	tempSource = zeros(leng,1);
	if(numOfSourcesGUI > 0)
		tempSource = conv(source_speech{1,1}, computed_rir.RIR_sources(:,i,1));
	end
	if(numOfSourcesGUI > 1)
		for j = 2:1:numOfSourcesGUI 
			tempSource = tempSource + conv(source_speech{1,1}, computed_rir.RIR_sources(:,i,j));
		end
	end 
	tempNoise = zeros(leng,1);
	if(numOfNoiseSourcesGUI > 0)
		tempNoise = conv(source_noise{1,1}, computed_rir.RIR_noise(:,i,1));
	end
	if(numOfNoiseSourcesGUI > 1)
		for k = 2:1:numOfNoiseSourcesGUI
			tempNoise = tempNoise + conv(source_noise{1,1}, computed_rir.RIR_noise(:,i,k));
		end
	end 
	
	mic(:,i) = tempSource + tempNoise;
	
	
	
end


fs = computed_rir.fs_RIR;
save('mic','mic','fs')

figure % to plot the two mic signals
hold on
plot(mic(:,1))
plot(mic(:,2))
hold off

 soundsc(mic(:,1),computed_rir.fs_RIR)



