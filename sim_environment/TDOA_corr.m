% TDOAgndTruth = sample delay between the direct path components of the two RIRs
computed_rir = load('Computed_RIRs.mat');

i = 1;
while (computed_rir.RIR_sources(i,1) <= 0)
	i = i+1;
end
first_index = i;
i = 1;
while (computed_rir.RIR_sources(i,2) <= 0)
	i = i+1;
end

TDOAgndTruth = abs(i-first_index);