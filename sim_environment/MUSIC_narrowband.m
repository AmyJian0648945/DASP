clear all;
computed_rir = load('Computed_RIRs.mat');

% Check if the sampling frequency is 44.1 kHz 
if(computed_rir.fs_RIR ~= 44.1e3)
	error('Error: sampling frequency is not 44.1kHz \n')
end 

% User defined variables (!!!)
lenMicSig = 10;		% length of desired microphone signal [sec]
mic = computeMicSig(computed_rir,lenMicSig);
numOfMics = size(computed_rir.RIR_sources,2);

%% STFT, with DFT window length L = 1024, and 50% overlap (each column contains stft estimate of each window) 
L = 1024;

for i=1:1:numOfMics
    stftMat(:,:,i) = spectrogram(mic(:,i), hamming(L), L./2, L);
end

