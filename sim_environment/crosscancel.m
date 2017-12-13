clear all;
load('Computed_RIRs.mat');
computed_rir = load('Computed_RIRs.mat'); 
load('HRTF.mat');
%% Variables Init
lenMicSig = 10; % 10 sec of speech
fs_resample = 8e3; % sampling frequency 
samplesToKeep = fs_resample.*lenMicSig;

RIR_length = 1600;
mic_synth = [1 4];
mic_listen = [1 4];

%% Speech 1 Reformatting
source_speech = cell(1,2); [source_speech{1,1},source_speech{1,2}] = audioread('speech1.wav');
source_speech{1,1} = resample(source_speech{1,1},fs_resample,source_speech{1,2}); %% resampling
source_speech1 = source_speech{1,1}(1:samplesToKeep); %% truncation

%%
RIR_sources = RIR_sources(1:RIR_length,:,:); %% truncating 1500 samples 
Lh = size(RIR_sources,1);
M = size(s_pos,1); %% number of loudspeakers
Lg = ceil(2*(Lh-1)./(M-2)); %% condition for having same amount of equations as of unknowns


temp = RIR_sources(:,mic_synth(2),1);
n = 3;
switch n
    case 1
        xL = zeros(Lg+Lh-1,1); xL(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xR = zeros(Lg+Lh-1,1); xR(1:RIR_length) = temp;
        
    case 2
        xL = zeros(Lg+Lh-1,1); xL(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xR = zeros(Lg+Lh-1,1); xR(1:RIR_length) = 0.2.*temp;
        
    case 3
        xL = zeros(Lg+Lh-1,1); xL(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xR = zeros(Lg+Lh-1,1); xR(1:RIR_length) = [zeros(10,1); temp(1:end-10)];
		
    case 4
        xL = zeros(Lg+Lh-1,1); xL(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xR = zeros(Lg+Lh-1,1); xR(1:RIR_length) = [zeros(10,1); temp(1:end-10)];
        xL = fftfilt(HRTF(:,1),xL,Lg+Lh-1);
        xR = fftfilt(HRTF(:,2),xR,Lg+Lh-1);
		
end


delta = ceil(sqrt(room_dim(1)^2 + room_dim(2)^2)*fs_RIR./340);

xL = [zeros(delta,1); xL(1:end-delta)];
xR = [zeros(delta,1); xR(1:end-delta)];

if(size(xL,1) >= Lh)
    %error('length of HRTF is larger than Lh');
end
x = [xL; xR];
H = zeros(2*(Lg+Lh-1),M*Lg);

for k = 1:1:M
    C = [RIR_sources(:,mic_synth(1),k); zeros(Lg-1,1)];
    R = [C(1) zeros(1,Lg-1)];
    H_k_left = toeplitz(C,R);
    
    C = [RIR_sources(:,mic_synth(2),k); zeros(Lg-1,1)];
    R = [C(1) zeros(1,Lg-1)];
    H_k_right = toeplitz(C,R);

    H(:,1+(k-1)*size(H_k_left,2):k*size(H_k_left,2)) = [H_k_left ; H_k_right];
end

%% Removing 0 rows
index_rows_zeros = (H==0);
sum_of_index = sum(index_rows_zeros,2);
index_2 = (sum_of_index==M*Lg);
H(index_2,:) =[]; 
x(index_2) = [];

%% Solving the SOE
g = H\x;
synth_error = norm(H*g-x)

%% Adding white noise to H
H = H + std(H(:,1)).* 0.05.* randn(size(H,1),size(H,2));

%% Filtering 

SourceFile = {'speech1.wav'};%, 'speech2.wav'};
NoiseFile = {}; %'White_noise1.wav', 'Babble_noise1.wav'};
flag_output = 1;
flag_input = 1;

mic = []; 
[mic, ~, ~] = computeMicSig4(computed_rir,10,flag_output,flag_input,SourceFile, NoiseFile,g,Lg, mic_listen); 

binaural_sig_expected = [fftfilt(xL,source_speech1,samplesToKeep) fftfilt(xR,source_speech1,samplesToKeep) ];
binaural_sig_with_g = [mic(:,1) mic(:,2)];

%% Sweet Spot calculation

radiusOfSweetSpot = 3e8 * M / (fs_resample*2*exp(1)*pi)

%% Plot results 
% figure 
% plot(x,'r');
% hold on 
% plot(H*g,'b');
% 
% figure, spy(H)
% soundsc(binaural_sig,8e3)









