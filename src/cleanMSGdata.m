%% cleanMSGdata
% helper function for webers_generally_observed

function cdata = cleanMSGdata(cdata)

% unpack
v2struct(cdata);

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);


AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% throw our bad traces

bad_trials = (sum(fA) == 0 | isnan(sum(fA)));
PID(:,bad_trials) = [];
LFP(:,bad_trials) = [];
fA(:,bad_trials) = [];
fB(:,bad_trials) = [];
paradigm(bad_trials) = [];
fly(bad_trials) = [];
orn(bad_trials) = [];

% extract filters and find gain for the firing rate
a = 35e3; z = 55e3;
try
	[K2,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);
catch
end

% filter the LFP
for i = 1:width(LFP)
	LFP(:,i) = LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,LFP(:,i));
	LFP(:,i) = LFP(:,i)*10; % to get the units right, now in mV
end

% extract filters and find gain for the LFP
[K1,LFP_pred,LFP_gain] = extractFilters(PID,LFP,'use_cache',true,'a',a,'z',z);

% repack everything
clear cdata
cdata = v2struct(a,z,K1,K2,PID,LFP,fA,fB,LFP_pred,fA_pred,fA_gain,LFP_gain,orn,fly,paradigm);
