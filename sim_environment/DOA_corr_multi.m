clear all; 
%%% DESCRIPTION:
% Cross-correlation-based direction-of-arrival (DOA) estimation for
% multiple sources.

computed_rir = load('Computed_RIRs.mat');

delay = TDOA_corr_separateSource();
numOfSources = size(computed_rir.RIR_sources,3);

d = norm(computed_rir.m_pos(1,:) - computed_rir.m_pos(2,:)); %intermicrophone distance
c = 340; %[m/s]

DOA_est = ones(1,numOfSources);
for i=1:1:numOfSources
	DOA_est(i) = delay(i) .* c ./ (computed_rir.fs_RIR .* d);
	if(DOA_est(i) < -1) DOA_est(i) = -1;%% to avoid complex numbers
	elseif(DOA_est(i) > 1) DOA_est(i) = 1;%% to avoid complex numbers
	else 
	end 
	DOA_est(i) = rad2deg(acos(DOA_est(i)));

end

save('DOA_est.mat','DOA_est')