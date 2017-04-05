% fig3_kinetics_whiff_based.m
% 
% created by Srinivas Gorur-Shandilya at 2:14 , 07 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;
dm = dataManager;

% main figure
main_fig = figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 4:-1:1
	ax(i) = subplot(2,2,i); hold on
end

% if you change this, you also have to change subplots
lag_vs_stim_plot = ax(3);
lag_vs_time_plot = ax(4);
gain_vs_mean_stim_plot = ax(1);
rho_vs_history_length_plot = ax(2);

% handle the plot within a plot for the gain vs. mu plot
subplots(1) = axes('Parent',main_fig);
subplots(2) = axes('Parent',main_fig);
subplots(3) = axes('Parent',main_fig);
subplots(4) = axes('Parent',main_fig);

subplots(1).Position = [.12 .75 .17 .17];
subplots(2).Position = [.3 .75 .17 .17];
subplots(3).Position = [.12 .56 .17 .17];
subplots(4).Position = [.3 .56 .17 .17];

gain_vs_mean_stim_plot.YColor = 'w';
gain_vs_mean_stim_plot.XColor = 'w';
gain_vs_mean_stim_plot.Position(1) = .08;
gain_vs_mean_stim_plot.Position(2) = .53;
gain_vs_mean_stim_plot.Position(3) = .4;
gain_vs_mean_stim_plot.XTick = [];
gain_vs_mean_stim_plot.YTick = [];



% supp figure
supp_fig = figure('outerposition',[0 0 1400 712],'PaperUnits','points','PaperSize',[1400 712]); hold on
for i = 10:-1:1
	axs(i) = subplot(2,5,i); hold on
end
delete(axs(5)); delete(axs(10));
axs(5) = []; axs(end) = [];





% define colors, etc. 
c = lines(10);
LFP_color = c(4,:);
firing_color = c(5,:);
model_color = [1 0 0];


% ##    ##    ###    ######## ##     ## ########     ###    ##       
% ###   ##   ## ##      ##    ##     ## ##     ##   ## ##   ##       
% ####  ##  ##   ##     ##    ##     ## ##     ##  ##   ##  ##       
% ## ## ## ##     ##    ##    ##     ## ########  ##     ## ##       
% ##  #### #########    ##    ##     ## ##   ##   ######### ##       
% ##   ### ##     ##    ##    ##     ## ##    ##  ##     ## ##       
% ##    ## ##     ##    ##     #######  ##     ## ##     ## ######## 

%  ######  ######## #### ##     ## ##     ## ##       #### 
% ##    ##    ##     ##  ###   ### ##     ## ##        ##  
% ##          ##     ##  #### #### ##     ## ##        ##  
%  ######     ##     ##  ## ### ## ##     ## ##        ##  
%       ##    ##     ##  ##     ## ##     ## ##        ##  
% ##    ##    ##     ##  ##     ## ##     ## ##        ##  
%  ######     ##    #### ##     ##  #######  ######## #### 


% analyse kinetics of LFP and firing rate during the naturalistic stimulus presentation
% load(getPath(dataManager,'aeb361c027b71938021c12a6a12a85cd'),'-mat');

load('/local-data/DA-paper/data-for-paper/nat-stim/v4_nat_stim.ORNData','-mat')
temp = od([1 3]);
clear od;
load('/local-data/DA-paper/data-for-paper/nat-stim/v3_nat_stim.ORNData','-mat')
od([1 6]) = [];
od = [temp od];
clear temp


% % first, fit the DA model to the data. (this was pre-fit, we're simply loading the fit here)
% load(getPath(dataManager,'b1b883899ab8c5ce1aed465819e75fce'));

% % generate DA model responses
% dd = ORNData;
% for i = 1:length(od)
% 	dd(i).stimulus = nanmean(od(i).stimulus,2);
% 	dd(i).firing_rate = DAModelv2(dd(i).stimulus,p(i));
% end

l_s = {'-','-','--','--','--','--','--','--','--','--'};
min_acceptable_corr = .5;
min_acceptable_lag = 5;
max_appeptable_lag = 300;
window_size = 1e3;
history_length = 1e3;
stim_thresh = .035;
clear l

for i = 1:length(od) % first one has stimulus which is too low, different paradigm 
	textbar(i,length(od))
	S = nanmean(od(i).stimulus,2); S = S - mean(S(1:5e3));
	R = nanmean(od(i).firing_rate,2);
	X = -nanmean(od(i).LFP,2);

	% find all whiffs
	[ons,offs] = findWhiffs(S);
	rm_this = offs + window_size > length(S) | ons - history_length < 1;
	ons(rm_this) = [];
	offs(rm_this) = [];


	lag_X = NaN*ons;
	lag_R = NaN*ons;
	time_since_last_whiff = NaN*ons;
	mean_s = NaN*ons;
	for j = 1:length(ons)
		s = S(ons(j):offs(j)+window_size); s = s - mean(s); s = s/std(s);
		x = X(ons(j):offs(j)+window_size); x = x - mean(x); x = x/std(x);
		r = R(ons(j):offs(j)+window_size); r = r - mean(r); r = r/std(r);
		lag_X(j) = finddelay(s,x);
		lag_R(j) = finddelay(s,r);
		mean_s(j) = mean(S(ons(j)-history_length:ons(j)));
	end


	rm_this = lag_X < min_acceptable_lag;
	% plot
	[~,data] = plotPieceWiseLinear(mean_s(~rm_this),lag_X(~rm_this),'nbins',19,'make_plot',false);
	plot(lag_vs_stim_plot,data.x,data.y,'Color',LFP_color,'LineStyle',l_s{i});

	% [~,data] = plotPieceWiseLinear(time_since_thresh_crossing,lag,'nbins',19,'make_plot',false);
	% plot(lag_vs_time_plot,data.x,data.y,'Color',firing_color,'LineStyle',l_s{i});

 % 	% LFP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	% [lag, ~, max_corr] = findLagAndMeanInWindow(S,X,window_size,25);
	% mean_x = circshift(vectorise(computeSmoothedStimulus(S,window_size)),history_length);
	% time_since_thresh_crossing = circshift(findTimeSinceThresholdCrossing(S,stim_thresh),-window_size);

	% % first strip out the NaNs
	% rm_this = isnan(lag);
	% lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];
	% time_since_thresh_crossing(rm_this) = [];

	% % then throw out some shitty data
	% rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr | lag > max_appeptable_lag;

	% [~,data] = plotPieceWiseLinear(mean_x(~rm_this),lag(~rm_this),'nbins',19,'make_plot',false);
	% plot(lag_vs_stim_plot,data.x,data.y,'LineStyle',l_s{i},'Color',LFP_color)

	% rm_this = lag < min_acceptable_lag | max_corr < min_acceptable_corr | time_since_thresh_crossing < 10 | lag > max_appeptable_lag;
	% [~,data] = plotPieceWiseLinear(time_since_thresh_crossing(~rm_this),lag(~rm_this),'make_plot',false,'nbins',19);
	% plot(lag_vs_time_plot,data.x,data.y,'LineStyle',l_s{i},'Color',LFP_color)

end

return


% labels -- main figure
xlabel(lag_vs_stim_plot,'\mu_{Stimulus} in preceding 300 ms (V)')
ylabel(lag_vs_stim_plot,'Lag (ms)')
set(lag_vs_stim_plot,'YLim',[0 160],'XLim',[1e-3 3],'XScale','log','YScale','linear','XTick',[1e-4 1e-3 1e-2 1e-1 1],'XMinorTick','off')


set(lag_vs_time_plot,'YLim',[0 160],'XLim',[10 1e4],'XScale','log')
xlabel(lag_vs_time_plot,'Time since odor encounter (ms)')
ylabel(lag_vs_time_plot,'Lag (ms)')


% labels -- supp. figure
xlabel(axs(1),'\mu_{Stimulus} in preceding 200ms (V)')
ylabel(axs(1),'Lag (ms)')
set(axs(1),'YLim',[0 140],'XLim',[0 0.6])

set(axs(5),'YLim',[0 140],'XLim',[10 5000],'XScale','log')
xlabel(axs(5),'Time since odor encounter (ms)')
ylabel(axs(5),'Lag (ms)')

% fake a fake legend for visual niceness
clear l
l(1) = plot(lag_vs_stim_plot,NaN,NaN,'Color',LFP_color);
l(2) = plot(lag_vs_stim_plot,NaN,NaN,'Color',firing_color);
L1 = legend(l,{'LFP','Firing rate'},'Location','southeast');

clear l
l(1) = plot(lag_vs_time_plot,NaN,NaN,'k-');
l(2) = plot(lag_vs_time_plot,NaN,NaN,'k--');
L2 = legend(l,{'ab2A','ab3A'},'Location','northwest');
L2.Position = [0.1700 0.4023 0.0800 0.0337];
uistack(lag_vs_time_plot,'top')


% ########    ###    ##     ##     ######      ###    #### ##    ## 
%    ##      ## ##   ##     ##    ##    ##    ## ##    ##  ###   ## 
%    ##     ##   ##  ##     ##    ##         ##   ##   ##  ####  ## 
%    ##    ##     ## ##     ##    ##   #### ##     ##  ##  ## ## ## 
%    ##    ######### ##     ##    ##    ##  #########  ##  ##  #### 
%    ##    ##     ## ##     ##    ##    ##  ##     ##  ##  ##   ### 
%    ##    ##     ##  #######      ######   ##     ## #### ##    ## 


% Now we compute the gain control timescale

history_lengths = round(logspace(1,log10(3e4),30)); % all the history lengths we look at, in ms
history_lengths(28) = [];
example_history_lengths = [3e2 1e3 3e3 1e4];
gain_mu = struct; gain_mu(length(example_history_lengths)).gain = []; gain_mu(1).mu = [];
rho = NaN(length(history_lengths),length(od));

use_this_segment = false(length(od(1).firing_rate),1);
use_this_segment(5e3:end-5e3) = true;


for i = 1:length(od)
	temp = od(i);
	pred = nanmean(temp.firing_projected,2); pred = pred(use_this_segment);
	resp = nanmean(temp.firing_rate,2);  resp = resp(use_this_segment);
	stim = nanmean(temp.stimulus,2); stim = stim - mean(stim(1:5e3));
	stim = stim(use_this_segment);

	% find when the valve opens
	[ons,offs] = findWhiffs(stim);
	ons = ons+finddelay(stim,resp);
	offs = offs+finddelay(stim,resp);

	% plot the gain in each of these windows
	[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

	rm_this = gain < 0 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];
	gain_err(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	if i > 3
		for j = 1:length(example_history_lengths)
			mu = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_lengths(j)));
			gain_mu(j).mu = [gain_mu(j).mu(:); mu(:)];
			gain_mu(j).gain = [gain_mu(j).gain(:); gain(:)];
		end
	end

	% also find rho for various values of the history length and plot it
	rho(:,i) = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

end

c = lines(10);
c(4:5,:) = [];
c(4,:) = [];
for i = 1:length(example_history_lengths)
	axes(subplots(i))
	plotPieceWiseLinear(gain_mu(i).mu,gain_mu(i).gain,'Color',c(i,:));
end


xlabel(gain_vs_mean_stim_plot,'\mu_{Stimulus} (V)','Color','k')
ylabel(gain_vs_mean_stim_plot,'ORN Gain (Hz/V)','Color','k')


% indicate these times on the fourth plot
for i = 1:length(example_history_lengths)
	plot(rho_vs_history_length_plot,[example_history_lengths(i) example_history_lengths(i)],[-1 ,1],'-','Color',c(i,:))
end

% first do ab2A
temp = errorbar(rho_vs_history_length_plot,history_lengths,nanmean(rho(:,1:2),2),nanstd(rho(:,1:2)/sqrt(2),[],2),'k');
temp = errorbar(rho_vs_history_length_plot,history_lengths,nanmean(rho(:,3:end),2),nanstd(rho(:,3:end)/sqrt(5),[],2),'k');
temp.LineStyle = '--';


set(rho_vs_history_length_plot,'XScale','log','YLim',[-1 1],'XTick',[10 1e2 1e3 1e4],'XLim',[10 3e4])
xlabel(rho_vs_history_length_plot,'Gain control timescale (ms)')
ylabel(rho_vs_history_length_plot,['Correlation between' char(10) 'gain and \mu_{stimulus}'])


prettyFig(main_fig,'fs',14);

% equalise axes
for i = 1:length(subplots)
	subplots(i).XLim = [2e-3 10];
	subplots(i).XTick = [ 1e-2 1e-1 1 1e1];
	subplots(i).YLim = [1e1 3e5];
	subplots(i).YTick = [10 100 1e3 1e4 1e5];
	subplots(i).Box = 'off';
	subplots(i).YScale = 'log';
	subplots(i).XScale = 'log';
	subplots(i).XMinorTick = 'off';
	subplots(i).YMinorTick = 'off';
end


subplots(2).YTickLabel = '';
subplots(2).XTickLabel = '';
subplots(1).XTickLabel = '';
subplots(4).YTickLabel = '';

L1.Box = 'off';
L2.Box = 'off';

if being_published
	snapnow
	delete(gcf)
end

return

%  ######  ##     ## ########  ########     ######## ####  ######   
% ##    ## ##     ## ##     ## ##     ##    ##        ##  ##    ##  
% ##       ##     ## ##     ## ##     ##    ##        ##  ##        
%  ######  ##     ## ########  ########     ######    ##  ##   #### 
%       ## ##     ## ##        ##           ##        ##  ##    ##  
% ##    ## ##     ## ##        ##           ##        ##  ##    ##  
%  ######   #######  ##        ##           ##       ####  ######   
   

% make space for the legend 
for i = [3 4 7 8]
	movePlot(axs(i),'right',.2)
end

% ##     ##  ######   ######   
% ###   ### ##    ## ##    ##  
% #### #### ##       ##        
% ## ### ##  ######  ##   #### 
% ##     ##       ## ##    ##  
% ##     ## ##    ## ##    ##  
% ##     ##  ######   ######   

% define what we want to work on
data_hashes = {'93ba5d68174e3df9f462a1fc48c581da','bcd4cf4fe12817d084a2b06f981161ee','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1','a33723c87a1216b274750734b4ee8820'};
odour_names = {'ethyl-acetate','1-pentanol','1-pentanol','2-butanone','isoamyl-acetate'};
orn_names = {'ab3A','ab3A','ab2A','ab2A','pb1A'};

markers = {'+','d','o','x','*'};
% core loop
for i = 1:length(data_hashes)-1
	clear cdata
	cdata = consolidateData2(getPath(dataManager,data_hashes{i}));
	if i < 4
		cdata.a = 25e3; cdata.z = 45e3;
	end
	cdata = cleanMSGdata(cdata,'extract_filter',false);

	plot_handles = plotMSGKinetics(cdata,axs(2));
	plot_handles(2).Color = firing_color;
	plot_handles(1).Color = LFP_color;

	% rescale the x axis
	temp = plot_handles(1).XData;
	temp = temp - min(temp);
	temp = temp/max(temp);
	plot_handles(1).XData =  temp;
	temp = plot_handles(2).XData;
	temp = temp - min(temp);
	temp = temp/max(temp);
	plot_handles(2).XData =  temp;

	plot_handles(1).MarkerSize = 10;
	plot_handles(2).MarkerSize = 10;

	plot_handles(1).Marker = markers{i};
	plot_handles(2).Marker = markers{i};
	t = [orn_names{i} char(10) odour_names{i}];
end

% fake some plots for a nice legend
clear l L
for i = 1:length(markers)-1
	l(i) = plot(axs(2),NaN,NaN,'Marker',markers{i},'Color','k','LineStyle','none');
	L{i} = [orn_names{i} ' ' odour_names{i}];
end
lh1 = legend(l,L,'Location','southeast');
lh1.Position = [0.43 0.7 0.12 0.1];
lh1.FontSize = 15;

set(axs(2),'YLim',[0 200],'XLim',[-0.1 1.1])
xlabel(axs(2),'Mean stimulus (rescaled)')

%  ######     ###    ########  ##        #######  ######## ########    ###    
% ##    ##   ## ##   ##     ## ##       ##     ##    ##       ##      ## ##   
% ##        ##   ##  ##     ## ##       ##     ##    ##       ##     ##   ##  
% ##       ##     ## ########  ##       ##     ##    ##       ##    ##     ## 
% ##       ######### ##   ##   ##       ##     ##    ##       ##    ######### 
% ##    ## ##     ## ##    ##  ##       ##     ##    ##       ##    ##     ## 
%  ######  ##     ## ##     ## ########  #######     ##       ##    ##     ## 

%% Now, show the timescale of gain control using Carlotta's data

%% global parameters
history_lengths = round(logspace(1.7,4,50)); % all the history lengths we look at, in ms
example_history_length = 300; % this history length shown in the first row, in ms

% load the data
if ~exist('orn_data','var')
	load(dm.getPath('86946ed05ec73186d8371166583141ba'))
end

do_these = [18 7 8 10 14 17 12];
odour_names = {'1-pentanol','methyl-butyrate','1-octen-3-ol','diethyl-succinate','ethyl-acetate','2-butanone','isoamyl-acetate'};

c = lines(length(do_these));
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

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(axs(6),history_lengths,rho,'-','Color',c(i,:));
end

set(axs(6),'XScale','log','YLim',[-1 0],'XTick',[1e1 1e2 1e3 1e4],'XLim',[10 1e4])
xlabel(axs(6),'Gain control timescale (ms)')
ylabel(axs(6),['Correlation between' char(10) 'gain and \mu_{stimulus}'])



% fake some plots for a nice legend
clear L l
for i = 1:length(do_these)
	l(i) = plot(axs(6),NaN,NaN,'Marker','o','MarkerFaceColor',c(i,:),'LineStyle','none');
	L{i} = [orn_data(do_these(i)).neuron_name ' ' odour_names{i}];
end
lh2 = legend(l,L,'Location','southeast');
lh2.Position = [0.43 0.15 0.15 0.2];
lh2.FontSize = 15;


% ########  ######## ########  ##       ####  ######     ###    ######## ########  ######  
% ##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##    ##       ##    ## 
% ##     ## ##       ##     ## ##        ##  ##        ##   ##     ##    ##       ##       
% ########  ######   ########  ##        ##  ##       ##     ##    ##    ######    ######  
% ##   ##   ##       ##        ##        ##  ##       #########    ##    ##             ## 
% ##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##    ##       ##    ## 
% ##     ## ######## ##        ######## ####  ######  ##     ##    ##    ########  ######  


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

	rm_this = gain<0.4 | gain_err < .8;
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

	plot(axs(7),history_lengths,rho,'.-','Color',c(i,:))
end

set(axs(7),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4])
set(axs(3),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.1 10])
xlabel(axs(3),['\mu_{Stimulus} in preceding ' oval(example_history_length) 'ms (norm)'])
ylabel(axs(3),'Gain (norm)')
xlabel(axs(7),'Gain control timescale (ms)')
ylabel(axs(7),['Correlation between' char(10) 'gain and \mu_{stimulus}'])
title(axs(3),['ab3A ORN,' char(10) ' methyl butyrate odorant'])

% fake some plots for a nice legend
clear l
for i = 1:length(c)
	l(i) = plot(axs(3),NaN,NaN,'Marker','o','MarkerFaceColor',c(i,:),'LineStyle','none','MarkerEdgeColor',c(i,:));
end
lh3 = legend(l,{'5/28','6/05','6/12','6/19','6/19'},'Location','southwest');



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

% tweak the model's tau_gain
tau_gain = [50 100 200 400];
stim = mean(orn_data(16).stimulus,2);
stim = stim(30e3:end);
stim = stim/nanmean(stim);


c = parula(length(tau_gain)+1);
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
	l(i) = plot(axs(4),mean_stim,gain(idx),'+-','Color',c(i,:));

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(axs(8),history_lengths,rho,'.-','Color',c(i,:))
end
set(axs(8),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4])
set(axs(4),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.05 10])
title(axs(4),['Simulations with ' char(10) 'varying \tau_{gain}'])
xlabel(axs(4),'\mu_{Stimulus} in preceding 300ms (V)')

clear l L
for i = 1:length(tau_gain)
	l(i) = plot(axs(4),NaN,NaN,'Marker','o','MarkerFaceColor',c(i,:),'LineStyle','none','MarkerEdgeColor',c(i,:));
	L{i} = ['\tau_{gain} = ' oval(4*tau_gain(i)) 'ms'];
end 
lh4 = legend(l,L);
lh4.Position = [0.855   0.62    0.035    0.08];
xlabel(axs(8),'Gain control timescale (ms)')

prettyFig(supp_fig,'fs',14);

if being_published
	snapnow
	delete(gcf)
end

% move some plots
axs(4).Position(1) = .82;
axs(8).Position(1) = .82;

axs(1).Position(1) = .12;
axs(5).Position(1) = .12;

%% Version Info
%
pFooter;




