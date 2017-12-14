clear all;
load('Computed_RIRs.mat');
computed_rir = load('Computed_RIRs.mat');
load('HRTF.mat');
%% Variables Init
lenMicSig = 10; % 10 sec of speech
fs_resample = 8e3; % sampling frequency
samplesToKeep = fs_resample.*lenMicSig;

RIR_length = 1600;
mic_synth = [1 2];
mic_listen = [1 2];

%% Speech 1 Reformatting
source_speech = cell(1,2); [source_speech{1,1},source_speech{1,2}] = audioread('speech1.wav');
source_speech{1,1} = resample(source_speech{1,1},fs_resample,source_speech{1,2}); %% resampling
source_speech1 = source_speech{1,1}(1:samplesToKeep); %% truncation

source_speechB = cell(1,2); [source_speechB{1,1},source_speechB{1,2}] = audioread('speech2.wav');
source_speechB{1,1} = resample(source_speechB{1,1},fs_resample,source_speechB{1,2}); %% resampling
source_speech1B = source_speechB{1,1}(1:samplesToKeep); %% truncation

%%
RIR_sources = RIR_sources(1:RIR_length,:,:); %% truncating 1500 samples
Lh = size(RIR_sources,1);
M = size(s_pos,1); %% number of loudspeakers
Lg = ceil(2*(Lh-1)./(M-2)); %% condition for having same amount of equations as of unknowns


temp = RIR_sources(:,mic_synth(2),1);
n = 7;
switch n
    case 1
        xL = zeros(Lg+Lh-1,1); xL(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xR = zeros(Lg+Lh-1,1); xR(1:RIR_length) = temp;
        
        xLB = zeros(Lg+Lh-1,1); xLB(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xRB = zeros(Lg+Lh-1,1); xRB(1:RIR_length) = temp;
        
    case 2
        xL = zeros(Lg+Lh-1,1); xL(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xR = zeros(Lg+Lh-1,1); xR(1:RIR_length) = 0.1.*temp;
        
        xLB = zeros(Lg+Lh-1,1); xLB(1:RIR_length) = 0.1.*RIR_sources(:,mic_synth(1),1);
        xRB = zeros(Lg+Lh-1,1); xRB(1:RIR_length) = RIR_sources(:,mic_synth(2),1);
        
    case 3
        xL = zeros(Lg+Lh-1,1); xL(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xR = zeros(Lg+Lh-1,1); xR(1:RIR_length) = [zeros(10,1); temp(1:end-10)];
        
        xLB = zeros(Lg+Lh-1,1); xLB(1:RIR_length) = [zeros(10,1); RIR_sources(1:end-10,mic_synth(1),1);];
        xRB = zeros(Lg+Lh-1,1); xRB(1:RIR_length) = RIR_sources(:,mic_synth(2),1);
        
    case 4
        xL = zeros(Lg+Lh-1,1); xL(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xR = zeros(Lg+Lh-1,1); xR(1:RIR_length) = RIR_sources(:,mic_synth(2),1);
        xL = fftfilt(HRTF(:,1),xL,Lg+Lh-1);
        xR = fftfilt(HRTF(:,2),xR,Lg+Lh-1);
        
        xLB = zeros(Lg+Lh-1,1); xLB(1:RIR_length) = RIR_sources(:,mic_synth(1),1);
        xRB = zeros(Lg+Lh-1,1); xRB(1:RIR_length) = RIR_sources(:,mic_synth(2),1);
        xLB = fftfilt(HRTF(:,2),xLB,Lg+Lh-1);
        xRB = fftfilt(HRTF(:,1),xRB,Lg+Lh-1);
    case 5
        xL = zeros(Lg+Lh-1,1); xL(1) = 1;
        xR = zeros(Lg+Lh-1,1); xR(1) = 1;
        
        xLB = zeros(Lg+Lh-1,1); xLB(1) =1;
        xRB = zeros(Lg+Lh-1,1); xRB(1) = 1;
    case 6
        xL = zeros(Lg+Lh-1,1); xL(1) = 1;
        xR = zeros(Lg+Lh-1,1); xR(1) = 0.5; %% attenuation ILD
        
        xLB = zeros(Lg+Lh-1,1); xLB(1) =0.5; %% attenuation ILD
        xRB = zeros(Lg+Lh-1,1); xRB(1) = 1;
    case 7
        xL = zeros(Lg+Lh-1,1); xL(1) = 1;
        xR = zeros(Lg+Lh-1,1); xR(4) = 1;%% delay of 3 samples ITD
        
        xLB = zeros(Lg+Lh-1,1); xLB(4) =1; %% delay of 3 samples ITD
        xRB = zeros(Lg+Lh-1,1); xRB(1) = 1;
    case 8
        Ltot =300;
        xL = zeros(Lg+Lh-1,1); xL(1:Ltot) = HRTF(1:Ltot,1);
        xR = zeros(Lg+Lh-1,1); xR(1:Ltot) = HRTF(1:Ltot,2);
        
        xLB = zeros(Lg+Lh-1,1); xLB(1:Ltot) =HRTF(1:Ltot,2); 
        xRB = zeros(Lg+Lh-1,1); xRB(1:Ltot) = HRTF(1:Ltot,1);
end


delta = ceil(sqrt(room_dim(1)^2 + room_dim(2)^2)*fs_RIR./340);

xL = [zeros(delta,1); xL(1:end-delta)];
xR = [zeros(delta,1); xR(1:end-delta)];
xLB = [zeros(delta,1); xLB(1:end-delta)];
xRB = [zeros(delta,1); xRB(1:end-delta)];

if(size(xL,1) >= Lh)
    %error('length of HRTF is larger than Lh');
end
x = [xL; xR];
H = zeros(2*(Lg+Lh-1),M*Lg);

xB = [xLB; xRB];

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
xB(index_2) = [];

%% Adding white noise to H
 H = H + std(H(:,1)).* 0.05.* randn(size(H,1),size(H,2));

%% Solving the SOE
g = H\x;
synth_error = norm(H*g-x)

gB = H\xB;
synth_errorB = norm(H*gB-xB)


%% Filtering

SourceFile = {'speech1.wav'};%, 'speech2.wav'};
NoiseFile = {}; %'White_noise1.wav', 'Babble_noise1.wav'};
flag_output = 1;
flag_input = 1;

mic = [];
[mic, ~, ~] = computeMicSig4(computed_rir,10,flag_output,flag_input,SourceFile, NoiseFile,g,Lg, mic_listen);

binaural_sig_expected = [fftfilt(xL,source_speech1,samplesToKeep) fftfilt(xR,source_speech1,samplesToKeep) ];
binaural_sig_with_g = [mic(:,1) mic(:,2)];

% FOR SPEECH 2
SourceFile = {'speech2.wav'};%, 'speech2.wav'};
NoiseFile = {}; %'White_noise1.wav', 'Babble_noise1.wav'};
flag_output = 1;
flag_input = 1;

micB = [];
[micB, ~, ~] = computeMicSig4(computed_rir,10,flag_output,flag_input,SourceFile, NoiseFile,gB,Lg, mic_listen);

binaural_sig_expectedB = [fftfilt(xLB,source_speech1B,samplesToKeep) fftfilt(xRB,source_speech1B,samplesToKeep) ];
binaural_sig_with_gB = [micB(:,1) micB(:,2)];

bin_expected = binaural_sig_expected + binaural_sig_expectedB;
bin_with_g = binaural_sig_with_g + binaural_sig_with_gB;

%% Sweet Spot calculation

radiusOfSweetSpot = 340 * M / (fs_resample*2*exp(1)*pi)

%% Plot results
figure
plot(x,'r');
hold on
plot(H*g,'b');
legend('x','H*g')
%
figure, spy(H)
%soundsc(binaural_sig,8e3)









