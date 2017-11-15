%%% DESCRIPTION: This function yields the psuedospectrum of a narrowband
%%%  signal
% 
%	Input: 
%		index_freq: the index of the frequency we're currently evaluating 
% 		freq:		the frequency bin that we're currently evaluating 
%		stftMat:	STFT matrix
%		computed_rir: microphone signals
%		theta:		angle data, i.e. 0:0.5:180
% 	Output: 
%		pw:			the pseudospectrum
% 



function [pw] = MUSIC_pseudoSpectrum_singFreq(index_freq,freq,stftMat,computed_rir,theta)

numOfMics = size(computed_rir.RIR_sources,2);
numOfSources = size(computed_rir.s_pos,1);



%% Create Ryy(w): MxM spatial correlation matrix

Ryy = zeros(numOfMics,numOfMics,size(stftMat, 2));

% Obtain 5x5 matrix for each frame (5 = # of mics)
for i=1:1:size(stftMat, 2) 
	
	% reshape to obtain a 5 x 1 vector instead of 1 x 1 x 5 vector
	y = reshape(stftMat(index_freq,i,:), [numOfMics,1]);
	
	% find autocorrelation matrix; Ryy: matrix of 5 x 5
	Ryy(:,:,i) = y*y';
end

% Take the average of all frames
Ryy = mean(Ryy, 3);



%% Create E(w): Mx(M-1) matrix

% Perform eig. decomposition
[E, eigVal] = eig(Ryy);
eigVal = diag(eigVal);
	
% Iteratively find the biggest # of eigenvalue and delete them, 
%  where # = numOfSources
for i=1:1:numOfSources
	[~,maxEigVal] = max(eigVal); % Find max eigenvalue
	eigVal(maxEigVal) = [];
	E(:,maxEigVal) = []; % Find and remove the corresponding eigenvector
end



%% Create g(w,theta) = array manifold vector; TDOA=tau, DOA=theta

dm = zeros(numOfMics,1);
c = 340; % [m/s]

% Calculate intermicrophone distance for all (compared to 1)
intermicDist = norm(computed_rir.m_pos(1,:) - computed_rir.m_pos(2,:)); 

for i=1:1:numOfMics
	dm(i) = (i-1).*intermicDist;
end

% Find tau(TDOA) function for each mic. Sampling frequency is removed from
%  the equation in order to yield units of [sec]
tau = -dm .* cos(theta) ./ c;

% Create g(w=max freq bin, theta)
w = freq .*2 .*pi;
g = exp(-1i .* w .* tau);

pw = diag(g'*E*E'*g).^-1; % should be scalar
end



