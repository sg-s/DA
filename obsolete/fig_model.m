
pHeader;

% uses dataManager for data integrity
dm = dataManager;

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on

%% grab MSG LFP data

opacity = .5;

[PID, LFP, fA, paradigm,~, ~, AllControlParadigms] = consolidateData(dm.getPath('bf79dfd769a97089e42beb0660174e84'),1);
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

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)) |  isnan(sum(LFP)));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];

% band pass all the LFP
try 
	load(dm.getPath('213e6122b7e0a414debcd5ded135ab20'),'filtered_LFP')
catch
	filtered_LFP = LFP;
	for i = 1:width(LFP)
		filtered_LFP(:,i) = 10*bandPass(LFP(:,i),1e4,Inf);
	end
end

% differentiate the LFP
dLFP = LFP;
for i = 1:length(paradigm)
	% first high pass them to remove spikes
	dLFP(:,i) = bandPass(LFP(:,i),Inf,10);
	dLFP(:,i) = -1e4*filtfilt(ones(10,1),10,[0; diff(dLFP(:,i))]);
end

% define limits on data
a = 10e3; z = 50e3;

% extract filters and compute gains
[K1,K1p,K1_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

% generate model responses
p.    r_b = 0.0134;
p.    r_d = 5.4150e-04;
p.theta_b = 0.0576;
p.theta_d = 1.4941;
p.      A = -40;

LFPm = LFP;
for i = 1:length(paradigm)
	LFPm(:,i) = simpleReceptorModelX(PID(:,i),p);
end

c = parula(max(paradigm)+1);
mean_stim = nanmean(PID(a:z,:)); 


% plot i/o curves for the model response. 
ax(2) = subplot(1,3,2); hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(LFPm(a:z,paradigm == i),2); y = y - nanmean(y);
	x = nanmean(K1p(a:z,paradigm == i),2);
	x = x - nanmean(x);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel(ax(2),'Projected Stimulus (V)')
ylabel(ax(2),'Model \DeltaLFP (mV)')

% show transduction gain
ax(3) = subplot(1,3,3); hold on
clear l
l(1) = plot(mean_stim,K1_gain,'k+');
l(2) = plot(ax(3),mean_stim,std(LFPm(a:z,:))./std(PID(a:z,:)),'r+');
set(ax(3),'XScale','log','YScale','log','XLim',[.1 3])
legend(l,{'Data','Model'})
xlabel(ax(3),'\mu_{Stimulus} (V)')
ylabel(ax(3),'Gain (mV/V)')

prettyFig
