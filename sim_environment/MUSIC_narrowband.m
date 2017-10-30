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

%% STFT, with DFT window length L = 1024, and 50% overlap
L = 1024;

% Make mic signal a multiple of 1024
% eucli_remainder = mod(length(mic), L);
% mic = [mic ; zeros((L - eucli_remainder),2)];
stftMatTemp = [];
stftMat = [];
for j=1:1:numOfMics
	i = 0;
	while(i.*L < length(mic))
		stftMatTemp = [stftMatTemp fft(mic(1+i.*L/2 : L + L./2.*i, j), L)];
		i = i+1;

	end 
	stftMat = cat(3, stftMat, stftMatTemp);
end


