clear all;

%% Running the create micSig script
SourceFile = {'speech1.wav'};%, 'speech2.wav'};
NoiseFile =  {'White_noise1.wav'};% {'White_noise1.wav','Babble_noise1.wav'};
computed_rir = load('Computed_RIRs.mat'); 
flag_output = 3;
flag_input = 4;
sourceLength = 3;

[mic, micSource, micNoise] = computeMicSig(computed_rir,sourceLength,flag_output,flag_input,SourceFile, NoiseFile); 
fs = computed_rir.fs_RIR;
load('SNR_in.mat');
%load('mic.mat');

%% WIDEBAND MUSIC Algorithm 
numOfMics = size(computed_rir.RIR_sources,2);
numOfSources = size(computed_rir.s_pos,1); %%% identifying Q 
theta = 0 : 0.5 : 180;
theta = deg2rad(theta);
%%% Find the STFT and the corresponding PSD 
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
%%% Computing P_hat(theta)
log_p_hat = 1./(L./2 -1) .* (sum(log(pw),2));
p_hat = exp((log_p_hat));
% find the location of all peaks
[peaks,locs] = findpeaks(abs(p_hat));
% sort the peaks in descending order, and only take into account the number
%  of peaks that correspond to the number of sources
[P,I] = sort(peaks,'descend');
DOA_est_rad = theta(locs(I(1:numOfSources)));
DOA_est = rad2deg(DOA_est_rad);
% Store the DOA estimate
save('DOA_est.mat','DOA_est');

% Plot the pseudospectrum p_hat
red_line = zeros(1,length(theta));
red_line(locs(I(1:numOfSources))) = abs(p_hat(locs(I(1:numOfSources))));
% Plot the pseudospectrum
figure('Name', 'Pseudospectrum \hat{p}(\theta)');
plot(abs(p_hat))
xlabel('Samples (0->180)')
ylabel('Amplitude')
title('Pseudospectrum using wideband MUSIC algorithm ')
hold on
stem(red_line)
%% DAS algorithm 
d = norm(computed_rir.m_pos(1,:) - computed_rir.m_pos(2,:));
c = 340; %% sound speed (m/s)
d_m = ((1:numOfMics)-1).*d;
delay = d_m*cos(DOA_est_rad)*fs./c %%% this are samples 
delay = round(delay)

%%%% STEP 1 : DELAY; cicrcular shifting for aligning signals 
for k=1:1:numOfMics
    mic(:,k) = circshift(mic(:,k),delay(k));
    micSource(:,k) = circshift(micSource(:,k),delay(k));
    micNoise(:,k) = circshift(micNoise(:,k),delay(k));
end

%%% STEP 2 : TRUNCATION; trucnating the max(delay) last/first samples 
maxDelay = max(delay)
if(maxDelay <= 0)
    micTrunc = mic(1:end+maxDelay,:);
else
    micTrunc = mic(maxDelay:end,:);
end

%%% STEP 3 : SUM; Adding signals %%%
DAS_out = sum(mic,2)./numOfMics;
speech_DAS = sum(micSource,2)./numOfMics;
noise_DAS = sum(micNoise,2)./numOfMics;

%% Plotting results
% figure 
% plot(mic(:,1))
% hold on 
% plot(DAS_out)

%% Computation of SNR for the DAS
%%% the theoretical gain should be M (#number of microphones) --> + 7dB for
%%% M = 5
VAD = abs(speech_DAS(:,1))>std(speech_DAS(:,1))*1e-3;
% Find the power of the signal of mic1 (where its active)
speechPower = var(speech_DAS(VAD==1,1));
noisePower = var(noise_DAS(:,1));
SNR_out_DAS = 10*log10(speechPower ./ noisePower);
save('SNR_out_DAS','SNR_out_DAS');
load('SNR_in.mat'); 

SNR
SNR_out_DAS
%% results 
%soundsc(DAS_out,fs)


