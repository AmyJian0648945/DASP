clear all;
%%% Group Members:
% Mattieu Fontaine r0692680
% Chin-Tze (Amy) Chien r0648945


computed_rir = load('Computed_RIRs.mat');

% Check if the sampling frequency is 44.1 kHz 
if(computed_rir.fs_RIR ~= 44.1e3)
	error('Error: sampling frequency is not 44.1kHz \n')
end 

% User defined variables (!!!)
lenMicSig = 1;		% length of desired microphone signal [sec]
mic = computeMicSig(computed_rir,lenMicSig);
numOfMics = size(computed_rir.RIR_sources,2);
numOfSources = size(computed_rir.s_pos,1);
theta = 0 : 0.5 : 180;
theta = deg2rad(theta);



%% Find the STFT and the corresponding PSD 
% DFT window length L = 1024, and 50% overlap (each column contains stft estimate of each window) 
%	stftMat:	matrix of L(frequency bins) x time bins x numOfMics 
%	freq:		matrix of L x 1
%	psd:		matrix of L x 1
L = 1024;
for i=1:1:numOfMics
    [stftMat(:,:,i), freq, ~, psd(:,:,i)] = spectrogram(mic(:,i), hann(L), [], 1*L, computed_rir.fs_RIR, 'psd');
end

% For each iterated frequency, evaluate the pseudo spectrum
%  where pw: matrix of theta x 1
pw = ones(length(theta),L./2);
for k=2:1:L./2
    freq_k = freq(k);
	pw(:,k-1) = MUSIC_pseudoSpectrum_singFreq(k,freq_k,stftMat,computed_rir,theta);
end



%% Computing P_hat(theta)

log_p_hat = 1./(L./2 -1) .* (sum(log(pw),2));
% for l=1:1:length(log_p_hat)
%     if (log_p_hat(l) > 75)
%       log_p_hat(l) = 75;
%     end
% end

p_hat = exp((log_p_hat));

% find the location of all peaks
[peaks,locs] = findpeaks(abs(p_hat));

% sort the peaks in descending order, and only take into account the number
%  of peaks that correspond to the number of sources
[P,I] = sort(peaks,'descend');
DOA_est = rad2deg(theta(locs(I(1:numOfSources))));

% Store the DOA estimate
save('DOA_est.mat','DOA_est');



%% Plotting results  

% Plot the pseudospectrum p_hat
figure('Name', 'Pseudospectrum \hat{p}(\theta)');
plot(theta, abs(p_hat))

% Plot the pseudospectrum 
figure('Name', 'Pseudospectrum p(\theta)');
plot(theta, abs(log(pw)))




