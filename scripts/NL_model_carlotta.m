
pHeader;

%% Can NL models explain what we see as fast gain control in Carlotta's data?
% First, a recap of the data: it is a binary stimulus, so should in principle be free of the problems that plagued the analysis of the naturalistic stimuli. When I performed the analysis of fast gain control on her data, I found a short timescale, and also found that the instatenous stimulus did not correlate with observed changed in gain:


%% global parameters
history_lengths = round(logspace(1.7,4,50)); % all the history lengths we look at, in ms
example_history_length = 300; % this history length shown in the first row, in ms

% load the data
if ~exist('orn_data','var')
	load(getPath(dataManager,'86946ed05ec73186d8371166583141ba'))
end

odour_names = {'1-pentanol','methyl-butyrate','1-octen-3-ol','diethyl-succinate','ethyl-acetate','2-butanone','isoamyl-acetate'};


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
axs(3) = subplot(1,2,1); hold on
axs(8) = subplot(1,2,2); hold on

do_these = [7 9 13 15 16];
c = lines(length(do_these));
clear l
for i = 1:length(do_these)
	temp = orn_data(do_these(i));
	pred = nanmean(temp.firing_projected,2); pred = pred(temp.use_this_segment);
	resp = nanmean(temp.firing_rate,2);  resp = resp(temp.use_this_segment);
	stim = nanmean(temp.stimulus,2); stim = stim(temp.use_this_segment);
 	stim = stim/nanmean(stim);

	% find when the valve opens
	[ons,offs] = findValveWhiffs(temp);

	% plot the gain in each of these windows
	[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

	gain = gain/nanmean(gain);

	rm_this = gain < 0 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	mean_stim = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_length));

	% % plot the gain vs. the mean stim after sorting it
	[mean_stim,idx] = sort(mean_stim);
	plot(axs(3),mean_stim,gain(idx),'+-','Color',c(i,:));

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(axs(8),history_lengths,rho,'.-','Color',c(i,:))
end

set(axs(3),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.1 10])
xlabel(axs(3),['\mu_{Stimulus} in preceding ' oval(example_history_length) ' ms (norm)'])
ylabel(axs(3),'ORN gain (norm)')

set(axs(8),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4],'XLim',[10 1.1e4])
xlabel(axs(8),'History length (ms)')

% fake some plots for a nice legend
clear l
for i = 1:length(c)
	l(i) = plot(axs(3),NaN,NaN,'Marker','o','MarkerFaceColor',c(i,:),'LineStyle','none','MarkerEdgeColor',c(i,:));
end
lh3 = legend(l,{'5/28','6/05','6/12','6/19','6/19'},'Location','southwest');

title(axs(3),['Binary stimulus' char(10) 'Experimental replicates'])


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now I fit a NLN model to this data. A NLN model fit to Carlotta's data also reproduces the curves which we think originate from fast gain control. (left column). How did I miss this earlier? It's really a function of the input nonlinearity -- a LN model does not show this effect. (right column), and that's what I used as my control. 

% generate synthetic data using the NLN model 
temp = orn_data(16);

clear data
resp = nanmean(temp.firing_rate,2);  resp = resp(temp.use_this_segment);
data.response = resp;
stim = nanmean(temp.stimulus,2); stim = stim(temp.use_this_segment);
stim = stim/nanmean(stim);
stim = stim - min(stim);
data.stimulus = [stim resp];

p.k_D = 0.3894;
p.n = 2.5;

R = NLNmodel(data.stimulus,p);

% now fit a linear model to this
K = fitFilter2Data(stim,R,'reg',1,'offset',200);
K = K(100:800);
filtertime = 1e-3*(1:length(K)) - .1;
time = 1e-3*(1:length(stim));
pred = convolve(time,stim,K,filtertime);

figure('outerposition',[0 0 800 801],'PaperUnits','points','PaperSize',[800 801]); hold on
axs(3) = subplot(2,2,1); hold on
axs(8) = subplot(2,2,3); hold on

resp = R;

% find when the valve opens
[ons,offs] = findValveWhiffs(temp);

% plot the gain in each of these windows
pred(pred<0) = 0;
[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

gain = gain/nanmean(gain);

rm_this = gain < 0 | gain_err < .8;
gain(rm_this) = [];
ons(rm_this) = [];
offs(rm_this) = [];

% find the mean stimulus in the preceding X ms in these windows
mean_stim = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_length));

% % plot the gain vs. the mean stim after sorting it
[mean_stim,idx] = sort(mean_stim);
plot(axs(3),mean_stim,gain(idx),'+-');
xlabel(axs(3),'Mean Stimulus')
ylabel(axs(3),'Gain')
title(axs(3),'NLN model')
set(axs(3),'XLim',[0.2 1.6],'YLim',[.3 1.6])

% also find rho for various values of the history length and plot it
rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

plot(axs(8),history_lengths,rho,'.-')
set(axs(8),'XScale','log','YLim',[-1 .4])
xlabel(axs(8),'History length (ms)')
ylabel(axs(8),'Spearman \rho')


% now repeat it, but use a LN model 
clear p
p.n = 2;
p.tau2 = 70;
p.tau1 = 20;
p.A = 0.3;
K0 = filter_gamma2(1:1e3,p);


R = filter(K0,1,stim); R(R<0) = 0;

% now fit a linear model to this
K = fitFilter2Data(stim,R,'reg',1,'offset',100);
K = K(50:800);
filtertime = 1e-3*(1:length(K)) - .05;
time = 1e-3*(1:length(stim));
pred = convolve(time,stim,K,filtertime);

axs(3) = subplot(2,2,2); hold on
axs(8) = subplot(2,2,4); hold on

resp = R;

% find when the valve opens
[ons,offs] = findValveWhiffs(temp);

% plot the gain in each of these windows
pred(pred<0) = 0;
[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

gain = gain/nanmean(gain);

rm_this = gain < 0 | gain_err < .8;
gain(rm_this) = [];
ons(rm_this) = [];
offs(rm_this) = [];

% find the mean stimulus in the preceding X ms in these windows
mean_stim = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_length));

% % plot the gain vs. the mean stim after sorting it
[mean_stim,idx] = sort(mean_stim);
plot(axs(3),mean_stim,gain(idx),'+-');
xlabel(axs(3),'Mean Stimulus')
ylabel(axs(3),'Gain')
title(axs(3),'LN model')
set(axs(3),'XLim',[0.2 1.6],'YLim',[.3 1.6])

% also find rho for various values of the history length and plot it
rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

plot(axs(8),history_lengths,rho,'.-')
set(axs(8),'XScale','log','YLim',[-1 .4])
xlabel(axs(8),'History length (ms)')
ylabel(axs(8),'Spearman \rho')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;


