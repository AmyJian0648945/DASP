clear all;
%%% DESCRIPTION:
%	Purpose of this file: to read input files, then truncate the file by the user specified lenghth, 
%	and resample the signal according to the sampling frequency of the RIR



computed_rir = load('Computed_RIRs.mat'); 

% Create the mic signal 
mic = computeMicSig(computed_rir); 

fs = computed_rir.fs_RIR;
save('mic','mic','fs')

figure % to plot the two mic signals
hold on
plot(mic(:,1))
plot(mic(:,2))
hold off

soundsc(mic(:,1),computed_rir.fs_RIR)



