clear all;
load('HRTF.mat')
%% Variables Init
lenMicSig = 10; % 10 sec of speech
fs_resample = 8e3; % sampling frequency 
samplesToKeep = fs_resample.*lenMicSig;
%% Speech 1
source_speech = cell(1,2); [source_speech{1,1},source_speech{1,2}] = audioread('speech1.wav');
source_speech{1,1} = resample(source_speech{1,1},fs_resample,source_speech{1,2}); %% resampling
x = source_speech{1,1}(1:samplesToKeep); %% truncation
binaural_sig1_1 = [x x];
binaural_sig2_1 = [x 0.5.*x];
right_sig = [0; 0; 0; x(1:end-3)];
binaural_sig3_1 = [x right_sig];
binaural_sig4_1 = [fftfilt(HRTF(:,1),x,samplesToKeep) fftfilt(HRTF(:,2),x,samplesToKeep) ];
%% Speech 2
source_speech2 = cell(1,2); [source_speech2{1,1},source_speech2{1,2}] = audioread('speech2.wav');
source_speech2{1,1} = resample(source_speech2{1,1},fs_resample,source_speech2{1,2}); %% resampling
x = source_speech2{1,1}(1:samplesToKeep); %% truncation
binaural_sig1_2 = [x x];
binaural_sig2_2 = [0.5.*x x];
left_sig = [0; 0; 0; x(1:end-3)];
binaural_sig3_2 = [left_sig  x];
binaural_sig4_2 = [fftfilt(HRTF(:,1),x,samplesToKeep) fftfilt(HRTF(:,2),x,samplesToKeep) ];