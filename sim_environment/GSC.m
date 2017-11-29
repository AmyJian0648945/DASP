clear all;

%% Create Speech Ref: Running DAS_BF file
DAS_BF % this loads all of the audio files, sets the flags...etc. 

%% Create Noise Ref: Blocking Matrix & NLMS

%%% Step 1: Blocking Matrix (Griffin-Jim)
% Create Blocking Matrix
Ca = [ones(numOfMics-1,1) diag(-ones(numOfMics-1,1))];
 
% flip micTrunc so that it matches the matrix on the text
micTrunc = micTrunc';

% y_1:M[k] = Ca * y_1:M[k] i.e. matrix multiplication: Ca*micTrunc(:,i)
% then concatenated. x_blockMatOuput: (M-1)xL
x_blockMatOuput = Ca * micTrunc;


%%% Step 2: NLMS, the multichannel adaptive filter (u=0.1, L=1024)








