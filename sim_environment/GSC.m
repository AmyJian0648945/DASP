clear all;

%% Create Speech Ref: Running DAS_BF file
DAS_BF % this loads all of the audio files, sets the flags...etc. 

%% Running GSC: NLMS & Updating W

% Create Blocking Matrix
Ca = [ones(numOfMics-1,1) diag(-ones(numOfMics-1,1))];

% Initialisations
u = 0.1; % step size
L = 1024; % filter size
X_blockL = zeros(numOfMics-1, L); % X_blockL: (M-1)xL
W_blockL = zeros(numOfMics-1, L); % W_blockL: (M-1)xL
micTrunc = micTrunc';
GSC_out = zeros(1,size(micTrunc,2));

listenToCa = [];
% Multichannel NLMS Function
for i=1:1:size(micTrunc, 2)
	
	% Check if there is speech signal present
	if(VAD==0)
		% Get windowed signal: MxL
		if(i < L) %i.e. if the window isn't within the signal, append with zeros
			X_blockL = [zeros(size(micTrunc,1), L-i) micTrunc(:,1:i)];
		else
			X_blockL = micTrunc(:, i-L+1:i);
		end 

		% Pass through blocking matrix
		X_blockL = Ca * X_blockL;
		listenToCa = [listenToCa X_blockL];

		% Calculate output: desiredSigal - filteredX
		GSC_out(i) = DAS_out(i) - sum(diag(X_blockL*W_blockL'));

		% Update Filter Coefficient (using previously calculated output): equ#5
		W_blockL = W_blockL + X_blockL.*(u .* GSC_out(i) ./  (norm(X_blockL,'fro').^2)  );
	else
		GSC_out(i) = DAS_out(i);
	end
	
end
 
listenToCa = fliplr(listenToCa);

%% Computation of SNR for the GSC
GSC_out = GSC_out';

% Find the power of the signal of mic1 (where its active)
noisePower = var(GSC_out(VAD==0,1));
speechPower = var(GSC_out(VAD==1,1)) - noisePower;

SNR_out_GSC = 10*log10(speechPower ./ noisePower)
save('SNR_out_GSC','SNR_out_GSC');






%% Plotting results
figure 
hold on 
plot(mic(:,1),'b')
plot(DAS_out,'r')
plot(GSC_out,'g')
% plot(speech_DAS,'k')
hold off

% Listening to the difference
% soundsc(speech_DAS,computed_rir.fs_RIR)
% soundsc(mic(:,1),computed_rir.fs_RIR)
% soundsc(DAS_out,computed_rir.fs_RIR)
% soundsc(GSC_out,computed_rir.fs_RIR)





