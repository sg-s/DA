% fig_gain_control_fast_broadly_observed
% makes a figure for the paper showing that gain control is fast and is broadly observed


pHeader;

%% global parameters
history_lengths = round(logspace(1.7,4,50)); % all the history lengths we look at, in ms
example_history_length = 300; % this history length shown in the first row, in ms

figure('outerposition',[0 0 801 800],'PaperUnits','points','PaperSize',[801 800]); hold on


%  ######     ###    ########  ##        #######  ######## ########    ###    
% ##    ##   ## ##   ##     ## ##       ##     ##    ##       ##      ## ##   
% ##        ##   ##  ##     ## ##       ##     ##    ##       ##     ##   ##  
% ##       ##     ## ########  ##       ##     ##    ##       ##    ##     ## 
% ##       ######### ##   ##   ##       ##     ##    ##       ##    ######### 
% ##    ## ##     ## ##    ##  ##       ##     ##    ##       ##    ##     ## 
%  ######  ##     ## ##     ## ########  #######     ##       ##    ##     ## 

% ########     ###    ########    ###    
% ##     ##   ## ##      ##      ## ##   
% ##     ##  ##   ##     ##     ##   ##  
% ##     ## ##     ##    ##    ##     ## 
% ##     ## #########    ##    ######### 
% ##     ## ##     ##    ##    ##     ## 
% ########  ##     ##    ##    ##     ## 

clearvars -except being_published history_lengths example_history_length

% load the data
if ~exist('orn_data','var')
	load('/local-data/DA-paper/fig4/Carlotta_Data.mat')
end



% ########  #### ######## ########         #######  ########   #######  ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##     ## ##       
% ##     ##  ##  ######   ######          ##     ## ##     ## ##     ## ########   ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ##     ## ##   ##         ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##     ## ##     ## ##    ##  ##    ## 
% ########  #### ##       ##       ###     #######  ########   #######  ##     ##  ######  
do_these = [18 7 8 10 14 17];

ax(1) = subplot(2,2,1); hold on
ax(2) = subplot(2,2,3); hold on

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

	rm_this = gain<0.4 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	mean_stim = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_length));

	% % plot the gain vs. the mean stim after sorting it
	[mean_stim,idx] = sort(mean_stim);
	l(i) = plot(ax(1),mean_stim,gain(idx),'+-','Color',c(i,:));

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(ax(2),history_lengths,rho,'.-','Color',c(i,:))
end
set(ax(2),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4])
set(ax(1),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.1 10])
xlabel(ax(1),['\mu_{Stimulus} in preceding ' oval(example_history_length) 'ms (norm)'])
ylabel(ax(1),'Gain (norm)')
xlabel(ax(2),'Gain control timescale (ms)')
ylabel(ax(2),'Spearman''s \rho')
title(ax(1),'ab3A ORN, Different odorants')

lh1 = legend(l,{'1-pentanol','methyl butyrate','1-octen-3-ol','diethyl succinate','ethyl acetate','2-butanone'});
lh1.Position = [0.1417    0.6   0.2010    0.1];
lh1.FontSize = 12;

% ########  #### ######## ########         #######  ########  ##    ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##     ## ###   ## ##    ## 
% ##     ##  ##  ##       ##              ##     ## ##     ## ####  ## ##       
% ##     ##  ##  ######   ######          ##     ## ########  ## ## ##  ######  
% ##     ##  ##  ##       ##              ##     ## ##   ##   ##  ####       ## 
% ##     ##  ##  ##       ##       ###    ##     ## ##    ##  ##   ### ##    ## 
% ########  #### ##       ##       ###     #######  ##     ## ##    ##  ######  


do_these = [12 17 25];
orn_data(25).use_these_trials = 1:orn_data(25).n_trials;
orn_data(25).use_this_segment = true(length(orn_data(25).firing_rate),1);


ax(1) = subplot(2,2,2); hold on
ax(2) = subplot(2,2,4); hold on

c = lines(length(do_these));
clear l
for i = 1:length(do_these)
	temp = orn_data(do_these(i));
	pred = nanmean(temp.firing_projected,2); pred = pred(temp.use_this_segment);
	resp = nanmean(temp.firing_rate,2);  resp = resp(temp.use_this_segment);
	stim = nanmean(temp.stimulus,2); stim = stim(temp.use_this_segment);
 	stim = stim/nanmean(stim);

 	if i < 3
		% find when the valve opens
		[ons,offs] = findValveWhiffs(temp);
	else
		% no valve signal, use excursions to define windows
		[ons,offs] = computeOnsOffs(resp>.1*max(resp));
	end

	% plot the gain in each of these windows
	[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

	gain = gain/nanmean(gain);

	rm_this = gain<0.4 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	mean_stim = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_length));

	% % plot the gain vs. the mean stim after sorting it
	[mean_stim,idx] = sort(mean_stim);
	l(i) = plot(ax(1),mean_stim,gain(idx),'+-','Color',c(i,:));

	% also find rho for various values of the history length and plot it
	if i < 3
		rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);
		plot(ax(2),history_lengths,rho,'.-','Color',c(i,:))
	else
		history_lengths2 = history_lengths;
		history_lengths2(history_lengths2<300) = [];
		rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths2);
		plot(ax(2),history_lengths2,rho,'.-','Color',c(i,:))
	end

	
end

set(ax(2),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4])
set(ax(1),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.1 10])
xlabel(ax(1),['\mu_{Stimulus} in preceding ' oval(example_history_length) 'ms (norm)' ])
ylabel(ax(1),'Gain (norm)')
xlabel(ax(2),'Gain control timescale (ms)')
ylabel(ax(2),'Spearman''s \rho')
title(ax(1),'Different OR types')

lh2 = legend(l,{'pb1A - isoamyl acetate','ab3A - 2-butanone','ab2A - ethyl acetate'});
lh2.Position = [0.61    0.6066    0.1923    0.0564];

prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end

%% Supplementary Figure
% The point of this supplementary figure is to validate the analysis that we do. 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on

% ########  ######## ########  ##       ####  ######     ###    ######## ########  ######  
% ##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##    ##       ##    ## 
% ##     ## ##       ##     ## ##        ##  ##        ##   ##     ##    ##       ##       
% ########  ######   ########  ##        ##  ##       ##     ##    ##    ######    ######  
% ##   ##   ##       ##        ##        ##  ##       #########    ##    ##             ## 
% ##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##    ##       ##    ## 
% ##     ## ######## ##        ######## ####  ######  ##     ##    ##    ########  ######  


do_these = [7 9 13 15 16];

ax(1) = subplot(2,2,1); hold on
ax(2) = subplot(2,2,3); hold on

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

	rm_this = gain<0.4 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	mean_stim = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_length));

	% % plot the gain vs. the mean stim after sorting it
	[mean_stim,idx] = sort(mean_stim);
	l(i) = plot(ax(1),mean_stim,gain(idx),'+-','Color',c(i,:));

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(ax(2),history_lengths,rho,'.-','Color',c(i,:))
end
set(ax(2),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4])
set(ax(1),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.1 10])
xlabel(ax(1),['\mu_{Stimulus} in preceding ' oval(example_history_length) 'ms (norm)'])
ylabel(ax(1),'Gain (norm)')
xlabel(ax(2),'Gain control timescale (ms)')
ylabel(ax(2),'Spearman''s \rho')
title(ax(1),['ab3A ORN,' char(10) ' methyl butyrate odorant'])
lh = legend(l,{'5/28','6/05','6/12','6/19','6/19'});
lh.Location = 'southwest';


%  ######  #### ##     ## ##     ## ##          ###    ######## ####  #######  ##    ##  ######  
% ##    ##  ##  ###   ### ##     ## ##         ## ##      ##     ##  ##     ## ###   ## ##    ## 
% ##        ##  #### #### ##     ## ##        ##   ##     ##     ##  ##     ## ####  ## ##       
%  ######   ##  ## ### ## ##     ## ##       ##     ##    ##     ##  ##     ## ## ## ##  ######  
%       ##  ##  ##     ## ##     ## ##       #########    ##     ##  ##     ## ##  ####       ## 
% ##    ##  ##  ##     ## ##     ## ##       ##     ##    ##     ##  ##     ## ##   ### ##    ## 
%  ######  #### ##     ##  #######  ######## ##     ##    ##    ####  #######  ##    ##  ######  

% first, fit a DA model to the methyl butyrate ab3A data
clear p
p.   s0 = -0.0011;
p.  n_z = 2;
p.tau_z = 50; % rounded off
p.  n_y = 2;
p.tau_y = 15.4531;
p.    C = 0.0126;
p.    A = 2.5013e+04;
p.    B = 384.5000;

% generate synthetic data responses and do the analyses on that
ax(1) = subplot(2,2,2); hold on
ax(2) = subplot(2,2,4); hold on

% tweak the model's tau_gain
tau_gain = [50 100 200];
stim = mean(orn_data(16).stimulus,2);
stim = stim(30e3:end);
stim = stim/nanmean(stim);


c = parula(length(tau_gain));
clear l
for i = 1:length(tau_gain)
	p.tau_z = tau_gain(i);

	% generate responses
	resp = DAModelv2(stim,p);

	% fit a filter to this
	K = fitFilter2Data(stim,resp,'filter_length',1e3,'offset',200);
	K = K(100:end-100);
	filtertime = (1:length(K)) - 100;
	pred = convolve(1:length(stim),stim,K,filtertime);


	% plot the gain in each of these windows
	[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

	gain = gain/nanmean(gain);

	rm_this = gain<0.4 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	mean_stim = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_length));

	% % plot the gain vs. the mean stim after sorting it
	[mean_stim,idx] = sort(mean_stim);
	l(i) = plot(ax(1),mean_stim,gain(idx),'+-','Color',c(i,:));

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(ax(2),history_lengths,rho,'.-','Color',c(i,:))
end
set(ax(2),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4])
set(ax(1),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.1 10])
title(ax(1),['Simulations with ' char(10) 'varying \tau_{gain}'])
L = {};
for i = 1:length(tau_gain)
	L{i} = ['\tau_{gain} = ' oval(4*tau_gain(i)) 'ms'];
end 
lh = legend(l,L);
lh.Position = [0.5835    0.6029    0.1738    0.1133];
xlabel(ax(2),'Gain control timescale (ms)')
ylabel(ax(2),'Spearman''s \rho')


prettyFig('fs',18)

%% Version Info
%
pFooter;



