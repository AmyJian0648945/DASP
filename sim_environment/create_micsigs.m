% Purpose of this file: to read input files, then truncate the file by the user specified lenghth, 
% and resample the signal according to the sampling frequency of the RIR

computed_rir = load('Computed_RIRs.mat'); % always call this before you start, because you need the sampling frequency of the RIR

% User defined variables
lenMicSig = 5; % length of desired microphone signal in seconds
fs_source = [0 0 0 0];

[source_speech1,fs_source(1)] = audioread('speech1.wav');
[source_speech2,fs_source(2)] = audioread('speech2.wav');
[source_noise1,fs_source(3)] = audioread('Babble_noise1.wav');
[source_noise2,fs_source(4)] = audioread('White_noise1.wav');

source_speech1 = resample(source_speech1,computed_rir.fs_RIR,fs_source(1));
source_speech2 = resample(source_speech2,computed_rir.fs_RIR,fs_source(2));
source_noise1 = resample(source_noise1,computed_rir.fs_RIR,fs_source(3));
source_noise1 = resample(source_noise1,computed_rir.fs_RIR,fs_source(4));

samplesToKeep = fs_RIR.*lenMicSig;
source_speech1 = source_speech1(1 : samplesToKeep);
source_speech2 = source_speech2(1 : samplesToKeep);
source_noise1 = source_noise1(1 : samplesToKeep);
source_noise2 = source_noise2(1 : samplesToKeep);

% Add RIR noise + source 1 + source 2
mic = RIR_noise;
mic(:,1) = mic(:,1) + RIR_sources(:,1,1) + RIR_sources(:,1,2);
mic(:,2) = mic(:,2) + RIR_sources(:,2,1) + RIR_sources(:,2,2);
mic(:,3) = mic(:,3) + RIR_sources(:,3,1) + RIR_sources(:,3,2);

save('mic','mic','fs_RIR')

figure % to plot the two mic signals
hold on
plot(mic(:,1))
plot(mic(:,2))
hold off



