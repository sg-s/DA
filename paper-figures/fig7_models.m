%% fig7_models.m
% 


dm = dataManager;
pHeader;

% plots to make
figure('outerposition',[0 0 901 902],'PaperUnits','points','PaperSize',[901 902]); hold on
clear ax
ax.da_model_vs_data = subplot(3,3,1:2); hold on
ax.da_model_vs_ln = subplot(3,3,3); hold on
ax.lag_constant_1 = subplot(3,3,4); hold on
ax.lag_constant_2 = subplot(3,3,5); hold on
ax.stimulus = subplot(3,3,7); hold on
ax.response = subplot(3,3,8); hold on
ax.weber_lfp = subplot(3,3,9); hold on

% define colors, etc. 
c = lines(10);
LFP_color = c(4,:);
firing_color = c(5,:);
model_color = c(3,:);


% analyse kinetics of LFP and firing rate during the naturalistic stimulus presentation
load(dm.getPath('aeb361c027b71938021c12a6a12a85cd'),'-mat');

% first, fit the DA model to the data. (this was pre-fit, we're simply loading the fit here)
load(dm.getPath('b1b883899ab8c5ce1aed465819e75fce'));

% generate DA model responses
dd = ORNData;
for i = 1:length(od)
	dd(i).stimulus = od(i).stimulus;
	for j = 1:od(i).n_trials
		dd(i).firing_rate(:,j) = DAModelv2(dd(i).stimulus(:,j),p(i));
	end
end

% fit LN models to this data
od = fitNL(od,'firing_rate');

time = 1e-3*(1:length(od(1).stimulus));
example_orn = 5;
plot(ax.da_model_vs_data,time,nanmean(od(example_orn).firing_rate,2),'Color',firing_color)
plot(ax.da_model_vs_data,time,nanmean(dd(example_orn).firing_rate,2),'Color',model_color)
set(ax.da_model_vs_data,'XLim',[0 70],'YLim',[0 150])
xlabel(ax.da_model_vs_data,'Time (s)')
ylabel(ax.da_model_vs_data,'Firing rate (Hz)')

% how well does this do, compared to a LN model?
r2_DA = NaN(length(od),1);
r2_LN = NaN(length(od),1);
for i = 1:length(od)
	resp = od(i).firing_rate;
	rm_this = sum(resp) == 0 | isnan(sum(resp));
	resp(:,rm_this) = [];
	pred = od(i).firing_projected;
	rm_this = sum(pred) == 0;
	pred(:,rm_this) = [];
	r2_LN(i) = rsquare(nanmean(resp,2),nanmean(pred,2));

	pred = dd(i).firing_rate;
	rm_this = sum(pred) == 0;
	pred(:,rm_this) = [];
	r2_DA(i) = rsquare(nanmean(resp,2),nanmean(pred,2));
end

plot(ax.da_model_vs_ln,[0 1],[0 1],'k--')
plot(ax.da_model_vs_ln,r2_LN,r2_DA,'k+')
ylabel(ax.da_model_vs_ln,'r^2 (DA Model)')
xlabel(ax.da_model_vs_ln,'r^2 (LN Model)')
set(ax.da_model_vs_ln,'XLim',[0 1],'YLim',[0 1])

% now show that the lags are constant with the DA model

min_acceptable_corr = .5;
min_acceptable_lag = 2;
clear l
for i = 1:length(od)
	S = nanmean(od(i).stimulus,2); S = S - mean(S(1:5e3));
	R = nanmean(od(i).firing_rate,2);
	X = -nanmean(od(i).LFP,2);
	DA_R = nanmean(dd(i).firing_rate,2);


	mean_x = vectorise(computeSmoothedStimulus(S,200));
	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,R,1e3,25);

	time_since_thresh_crossing = findTimeSinceThresholdCrossing(S,mean(S));

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% then throw out some shitty data
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr | time_since_thresh_crossing < 10 | lag > 300;
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% plot
	axes(ax.lag_constant_1)
	l(2) = plotPieceWiseLinear(mean_x,lag,'Color',firing_color,'nbins',19);

	axes(ax.lag_constant_2)
	plotPieceWiseLinear(time_since_thresh_crossing,lag,'Color',firing_color,'nbins',19);

	% now do the same with the DA model prediction
	mean_x = vectorise(computeSmoothedStimulus(S,200));
	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,DA_R,1e3,25);

	time_since_thresh_crossing = findTimeSinceThresholdCrossing(S,mean(S));

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% then throw out some shitty data
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr | time_since_thresh_crossing < 10 | lag > 300;
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	axes(ax.lag_constant_1)
	l(3) = plotPieceWiseLinear(mean_x,lag,'Color',model_color,'nbins',19);

	axes(ax.lag_constant_2)
	plotPieceWiseLinear(time_since_thresh_crossing,lag,'Color',model_color,'nbins',19);

	mean_x = vectorise(computeSmoothedStimulus(S,200));
	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,X,1e3,25);
	time_since_thresh_crossing = findTimeSinceThresholdCrossing(S,mean(S));

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% then throw out some shitty data
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr | time_since_thresh_crossing < 10 | lag > 300;
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	time_since_thresh_crossing(rm_this) = [];

	% plot
	axes(ax.lag_constant_1)
	l(1) = plotPieceWiseLinear(mean_x,lag,'Color',LFP_color,'nbins',19);

	axes(ax.lag_constant_2)
	plotPieceWiseLinear(time_since_thresh_crossing,lag,'Color',LFP_color,'nbins',19);

end

% labels -- main figure
xlabel(ax.lag_constant_1,'\mu_{Stimulus} in preceding 200ms (V)')
ylabel(ax.lag_constant_1,'Lag (ms)')
set(ax.lag_constant_1,'YLim',[0 140],'XLim',[0 0.45])
L = legend(l,{'LFP','Firing Rate','DA Model'},'Location','southeast');
L.Position(1) = .7;

set(ax.lag_constant_2,'YLim',[0 140],'XLim',[10 1e4],'XScale','log','XTick',[10 1e2 1e3 1e4])
xlabel(ax.lag_constant_2,'Time since odor encounter (ms)')
ylabel(ax.lag_constant_2,'Lag (ms)')

% now the LFP model

% core parameters 
T = 20e3; % total length
pulse_on = 15e3;
pulse_off = 16e3;

clear p
p.A = 1e3;
p.B = 1;
p.tau = 1e-3;
p.ko = 1;


background_levels = logspace(-2,2,10);
c = parula(length(background_levels)+1);

% make sure we don't get negative responses
p.ko = -log(background_levels(1))*1.1;

all_gain = NaN*background_levels;
all_mu = NaN*background_levels;
for i = 1:length(background_levels)
	S = background_levels(i) + zeros(T,1);
	S(pulse_on:pulse_off) = 2*background_levels(i);
	plot(ax.stimulus,S,'Color',c(i,:));

	R = LFPmodel(S,p);
	% plot(ax(2),R,'Color',c(i,:))

	temp = R - mean(R(pulse_on-1e3:pulse_on-1));

	plot(ax.response,temp/max(temp),'Color',c(i,:))


	% compute gain 
	g = max(temp)/(max(S) - min(S));
	plot(ax.weber_lfp,min(S),g,'+','Color',c(i,:));
	all_gain(i) = g;
	all_mu(i) = max(S) - min(S);
end


set(ax.stimulus,'XLim',[length(S) - 10e3 length(S)])
set(ax.response,'XLim',[length(S) - 10e3 length(S)])

% draw a weber's fit to this
x = all_mu;
y = all_gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax.weber_lfp,sort(x),cf(sort(x)),'r')

set(ax.stimulus,'YScale','log')
ylabel(ax.stimulus,'Stimulus')
ylabel(ax.response,'r(t) (rescaled)')
set(ax.response,'YLim',[-.1 1])
set(ax.weber_lfp,'XScale','log','YScale','log','XTick',logspace(-2,2,5),'YTick',logspace(-2,2,5))
xlabel(ax.weber_lfp,'Background S')
ylabel(ax.weber_lfp,'Gain (\DeltaR / \DeltaS)')



prettyFig('fs',12)

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


