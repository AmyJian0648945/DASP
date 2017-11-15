clear all;
%%% DESCRIPTION:
%	Purpose of this file: to read input files, then truncate the file by the user specified lenghth, 
%	and resample the signal according to the sampling frequency of the RIR



computed_rir = load('Computed_RIRs.mat'); 

% Function call - obtain the mic signal
mic = computeMicSig(computed_rir,5); 

% Save mic signal as a .mat variable
fs = computed_rir.fs_RIR;
save('mic','mic','fs')

figure 
hold on
plot(mic(:,1))
plot(mic(:,2))
hold off

soundsc(mic(:,1),computed_rir.fs_RIR)



