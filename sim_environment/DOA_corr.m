clear all; 
%%% DESCRIPTION:
% Cross-correlation-based direction-of-arrival (DOA) estimation

computed_rir = load('Computed_RIRs.mat');
d = norm(computed_rir.m_pos(1,:) - computed_rir.m_pos(2,:)); %intermicrophone distance
c = 340; %[m/s]

delay = TDOA_corr();
DOA_est = delay .* c ./ (computed_rir.fs_RIR .* d);
if(DOA_est < -1) DOA_est = -1;
elseif(DOA_est > 1) DOA_est = 1;
else 
end 
DOA_est = rad2deg(acos(DOA_est));




save('DOA_est.mat','DOA_est')