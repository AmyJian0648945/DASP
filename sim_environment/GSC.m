clear all;

%% Create Speech Ref: Running DAS_BF file
DAS_BF % this loads all of the audio files, sets the flags...etc. 

%% Create Noise Ref: Blocking Matrix & NLMS

%%% Step 1: Blocking Matrix (Griffin-Jim)
% Create Blocking Matrix
Ca = [ones(numOfMics-1,1) diag(-ones(numOfMics-1,1))];

% y_1:M[k] = Ca * y_1:M[k] i.e. matrix multiplication: Ca*micTrunc(:,i)
% then concatenated. 
%	X: (M-1)xN
%	Ca: (M-1)xM
%	micTrunc: MxN
micTrunc = micTrunc'; % micTrunc: MxN
X = Ca * micTrunc;



%%% Step 2: NLMS, the multichannel adaptive filter (u=0.1, L=1024)
u = 0.1; % step size
L = 1024; % filter size
X = X'; % X: Nx(M-1), where N>>L
X_blockL = zeros(L,size(X,2)); % X_blockL: Lx(M-1)
W_blockL = zeros(L,size(X,2)); % W_blockL: Lx(M-1)
GSC_out = [];

if(mod(size(X,1),L) ~= 0)
	X = [X; zeros(L-mod(size(X,1),L), size(X,2))];
end

% delay the DAS BF
% DAS_out = [zeros(L./2,1) ; DAS_out]; 

% Multichannel NLMS update function
for i=1:L:size(X,1) %% Toggles filter block

	% Isolate the block for this filter
	X_blockL = X(i:i+L-1,:);
	if(i+L > size(DAS_out,1))
		
	else
		DAS_blockL = DAS_out(i:i+L-1,:);
	end
	
	
	% for each L block, we compute W:Lx(M-1). See equation (5) on exercise doc
	for j=1:1:L
		
		% zero out W_blockL?? 
		W_blockL = zeros(size(W_blockL,1), size(W_blockL,2));
		
		% Within square bracket
		tempSum = 0;
		for k=1:1:size(X,2)
			tempSum = tempSum + (W_blockL(:,k)') * X_blockL(:,k);
		end
		tempSum = DAS_blockL(j,:) - tempSum;
		
		% Outside of square bracket
		tempSum = (tempSum .* u ./ norm(X_blockL(j,:),'fro').^2 );
		tempUpdate = X_blockL(j,:) .* tempSum;
		if(j<L)
			W_blockL(j+1,:) = W_blockL(j,:) + tempUpdate;
		end 
	end
	
	% Speech ref - Noise ref 
	GSC_out = [GSC_out ; DAS_blockL-sum(W_blockL,2)];
	
	
end 

GSC_out = GSC_out(1:size(DAS_out,1),:);

%% Plotting results
close all;
figure 
hold on 
plot(mic(:,1),'b')
plot(DAS_out,'r')
plot(GSC_out,'g')
plot(speech_DAS,'k')
hold off

% Listening to the difference
% soundsc(mic(:,1),computed_rir.fs_RIR)
% soundsc(DAS_out,computed_rir.fs_RIR)
% soundsc(GSC_out,computed_rir.fs_RIR)
% soundsc(speech_DAS,computed_rir.fs_RIR)




