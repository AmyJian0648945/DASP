function [pw] = MUSIC_pseudoSpectrum_singFreq(index_freq)

% Check if the sampling frequency is 44.1 kHz 
if(computed_rir.fs_RIR ~= 44.1e3)
	error('Error: sampling frequency is not 44.1kHz \n')
end 

% User defined variables (!!!)
lenMicSig = 1;		% length of desired microphone signal [sec]
mic = computeMicSig(computed_rir,lenMicSig);
numOfMics = size(computed_rir.RIR_sources,2);
numOfSources = size(computed_rir.s_pos,1);



%% Find the STFT and the corresponding PSD 
% DFT window length L = 1024, and 50% overlap (each column contains stft estimate of each window) 
L = 1024;
for i=1:1:numOfMics
    [stftMat(:,:,i), freq, ~, psd(:,:,i)] = spectrogram(mic(:,i), L, [], 2.*L, computed_rir.fs_RIR, 'psd');
end

% Averaged over time and mic, find the max power; index is the frequency
% bin with the highest power
psd = mean(mean(psd, 2), 3); % avg over diff. frames & mic
[~, indexMaxFreq] = max(psd);
wMax = freq(index_freq) .* 2 .* pi;



%% Create Ryy(w): MxM spatial correlation matrix

Ryy = zeros(numOfMics,numOfMics,size(stftMat, 2));

% Obtain 5x5 matrix for each frame (5 = # of mics)
for i=1:1:size(stftMat, 2) 
	y = reshape(stftMat(index_freq,i,:), [numOfMics,1]);
	Ryy(:,:,i) = y*y';
end

% Take the average of all frames
Ryy = mean(Ryy, 3);



%% Create E(w): Mx(M-1) matrix

% Perform eig. decomposition, and find the largest eigenvalue
[E, eigVal] = eig(Ryy);
eigVal = diag(eigVal);
% [~,maxEigVal] = max(eigVal); % Find max eigenvalue
% E(:,maxEigVal) = []; % Find and remove the corresponding eigenvector
	
for i=1:1:numOfSources
	[~,maxEigVal] = max(eigVal); % Find max eigenvalue
	eigVal(maxEigVal) = [];
	E(:,maxEigVal) = []; % Find and remove the corresponding eigenvector
end


%% Create g(w,theta) = array manifold vector; TDOA=tau, DOA=theta

dm = zeros(numOfMics,1);
theta = 0 : 0.5 : 180;
theta = deg2rad(theta);
c = 340; % [m/s]

% Calculate intermicrophone distance for all (compared to 1)
intermicDist = norm(computed_rir.m_pos(1,:) - computed_rir.m_pos(2,:)); 

for i=1:1:numOfMics
	dm(i) = (i-1).*intermicDist;
end

% Find tau(TDOA) function for each mic
tau = -dm .* cos(theta) ./ c;

% Create g(w=max freq bin, theta)
g = exp(-1i .* wMax .* tau);



%% Evaluate the pseudospectrum

% Find the true DOA
pw = diag(g'*E*E'*g).^-1; % should be scalar
end



