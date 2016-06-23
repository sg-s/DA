% fig_kinetics.m
% 
% created by Srinivas Gorur-Shandilya at 2:14 , 07 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;
dm = dataManager;

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 4:-1:1
	ax(i) = subplot(2,2,i); hold on
end


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
load(dm.getPath('aeb361c027b71938021c12a6a12a85cd'),'-mat');

c = lines(3);

min_acceptable_corr = .5;
min_acceptable_lag = 2;
clear l
for i = 1:length(od)
	S = nanmean(od(i).stimulus,2); S = S - mean(S(1:5e3));
	R = nanmean(od(i).firing_rate,2);
	X = -nanmean(od(i).LFP,2);

	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,R,1e3,25);
	time_since_thresh_crossing = findTimeSinceThresholdCrossing(S,mean(S));
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
	lag(rm_this) = [];
	mean_x(rm_this) = [];
	t = time_since_thresh_crossing;
	t(rm_this) = []; t(t<10) = NaN;
	lag(lag>300) = NaN; % remove some obvious outliers
	axes(ax(1))
	l(2) = plotPieceWiseLinear(mean_x,lag,'Color',c(1,:),'nbins',19);

	axes(ax(3))
	plotPieceWiseLinear(t,lag,'Color',c(1,:),'nbins',19);


	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,X,1e3,25);
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
	lag(rm_this) = [];
	mean_x(rm_this) = [];
	t = time_since_thresh_crossing;
	t(rm_this) = []; t(t<10) = NaN;

	axes(ax(1))
	l(1) = plotPieceWiseLinear(mean_x,lag,'Color',c(2,:),'nbins',19);

	axes(ax(3))
	plotPieceWiseLinear(t,lag,'Color',c(2,:),'nbins',19);

end
xlabel(ax(1),'\mu_{Stimulus} in preceding 1s (V)')
ylabel(ax(1),'Lag (ms)')
set(ax(1),'YLim',[0 140],'XLim',[0 0.6])
L = legend(l,{'LFP','Firing Rate'},'Location','southeast');
title(ax(1),['ab3A' char(10) 'ethyl-acetate'])

set(ax(3),'YLim',[0 140],'XLim',[10 5000],'XScale','log')
xlabel(ax(3),'Time since odor encounter (ms)')
ylabel(ax(3),'Lag (ms)')


% Now we 

history_lengths = round(logspace(1.7,4,50)); % all the history lengths we look at, in ms
example_history_lengths = [100 10e3];
gain_mu = struct; gain_mu(length(example_history_lengths)).gain = []; gain_mu(1).mu = [];
rho = NaN(length(history_lengths),length(od));

use_this_segment = false(length(od(1).firing_rate),1);
use_this_segment(5e3:end-5e3) = true;

for i = [2 3 5 6]
	temp = od(i);
	pred = nanmean(temp.firing_projected,2); pred = pred(use_this_segment);
	resp = nanmean(temp.firing_rate,2);  resp = resp(use_this_segment);
	stim = nanmean(temp.stimulus,2); stim = stim - mean(stim(1:5e3));
	stim = stim(use_this_segment);

	% find when the valve opens
	[ons,offs] = findWhiffs(stim);

	% plot the gain in each of these windows
	[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

	rm_this = gain < 0 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	for j = 1:length(example_history_lengths)
		mu = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_lengths(j)));
		gain_mu(j).mu = [gain_mu(j).mu(:); mu(:)];
		gain_mu(j).gain = [gain_mu(j).gain(:); gain(:)];
	end

	% also find rho for various values of the history length and plot it
	rho(:,i) = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

end

errorbar(ax(4),history_lengths,nanmean(rho,2),nanstd(rho,[],2),'k')
set(ax(4),'XScale','log','YLim',[-1 .4],'XTick',[10 1e2 1e3 1e4],'XLim',[10 1e4])
xlabel(ax(4),'Gain control timescale (ms)')
ylabel(ax(4),['Correlation between' char(10) 'sgain and \mu_{stimulus}'])

% show the gain vs. the mean stimulus
c = lines(4);
axes(ax(2))
plotPieceWiseLinear(gain_mu(end).mu,gain_mu(end).gain,'nbins',10,'Color',c(3,:));

axes(ax(2))
plotPieceWiseLinear(gain_mu(1).mu,gain_mu(1).gain,'nbins',10,'Color',c(4,:));

% fit Weber's Law to this
x = gain_mu(1).mu;
y = gain_mu(1).gain;
options = fitoptions(fittype('poly1'));
options.Lower = [-1 -Inf];
options.Upper = [-1 Inf];
cf_temp = fit(log(x(:)),log(y(:)),'poly1',options);
cf = fit(x(:),y(:),'power1');
warning off
cf.a = exp(cf_temp.p2); cf.b = -1;
warning on
plot(ax(2),sort(x),cf(sort(x)),'r')



set(ax(2),'XScale','log','YScale','log')
xlabel(ax(2),'\mu_{Stimulus} (V)')
ylabel(ax(2),'Gain (Hz/V)')

prettyFig('fs',14);

if being_published
	snapnow
	delete(gcf)
end

%  ######  ##     ## ########  ########     ######## ####  ######   
% ##    ## ##     ## ##     ## ##     ##    ##        ##  ##    ##  
% ##       ##     ## ##     ## ##     ##    ##        ##  ##        
%  ######  ##     ## ########  ########     ######    ##  ##   #### 
%       ## ##     ## ##        ##           ##        ##  ##    ##  
% ##    ## ##     ## ##        ##           ##        ##  ##    ##  
%  ######   #######  ##        ##           ##       ####  ######   
    
figure('outerposition',[0 0 1200 712],'PaperUnits','points','PaperSize',[1200 712]); hold on
clear ax
for i = 8:-1:1
	ax(i) = subplot(2,4,i); hold on
end
delete(ax(4)); delete(ax(8))

% make space for the legend 
for i = [2 3 6 7]
	movePlot(ax(i),'right',.2)
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
	cdata = consolidateData2(dm.getPath(data_hashes{i}));
	if i < 4
		cdata.a = 25e3; cdata.z = 45e3;
	end
	cdata = cleanMSGdata(cdata,'extract_filter',false);

	plot_handles = plotMSGKinetics(cdata,ax(1));

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
	l(i) = plot(ax(1),NaN,NaN,'Marker',markers{i},'Color','k','LineStyle','none');
	L{i} = [orn_names{i} ' ' odour_names{i}];
end
lh1 = legend(l,L,'Location','southeast');
lh1.Position = [0.31 0.7 0.15 0.1];
lh1.FontSize = 15;

set(ax(1),'YLim',[0 200],'XLim',[-0.1 1.1])
xlabel(ax(1),'Mean stimulus (rescaled)')

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

	plot(ax(5),history_lengths,rho,'-','Color',c(i,:));
end

set(ax(5),'XScale','log','YLim',[-1 0],'XTick',[1e2 1e3 1e4],'XLim',[100 1e4])
xlabel(ax(5),'Gain control timescale (ms)')
ylabel(ax(5),['Correlation between' char(10) 'gain and \mu_{stimulus}'])



% fake some plots for a nice legend
clear L l
for i = 1:length(do_these)
	l(i) = plot(ax(5),NaN,NaN,'Marker','o','MarkerFaceColor',c(i,:),'LineStyle','none');
	L{i} = [orn_data(do_these(i)).neuron_name ' ' odour_names{i}];
end
lh2 = legend(l,L,'Location','southeast');
lh2.Position = [0.31 0.15 0.15 0.2];
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
	plot(ax(2),mean_stim,gain(idx),'+-','Color',c(i,:));

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(ax(6),history_lengths,rho,'.-','Color',c(i,:))
end

set(ax(6),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4])
set(ax(2),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.1 10])
xlabel(ax(2),['\mu_{Stimulus} in preceding ' oval(example_history_length) 'ms (norm)'])
ylabel(ax(2),'Gain (norm)')
xlabel(ax(6),'Gain control timescale (ms)')
ylabel(ax(6),['Correlation between' char(10) 'gain and \mu_{stimulus}'])
title(ax(2),['ab3A ORN,' char(10) ' methyl butyrate odorant'])

% fake some plots for a nice legend
clear l
for i = 1:length(c)
	l(i) = plot(ax(2),NaN,NaN,'Marker','o','MarkerFaceColor',c(i,:),'LineStyle','none','MarkerEdgeColor',c(i,:));
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
	l(i) = plot(ax(3),mean_stim,gain(idx),'+-','Color',c(i,:));

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(ax(7),history_lengths,rho,'.-','Color',c(i,:))
end
set(ax(7),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4])
set(ax(3),'XScale','log','YScale','log','XLim',[.1 2],'YLim',[.05 10])
title(ax(3),['Simulations with ' char(10) 'varying \tau_{gain}'])
xlabel(ax(3),'\mu_{Stimulus} in preceding 300ms (V)')

clear l L
for i = 1:length(tau_gain)
	l(i) = plot(ax(3),NaN,NaN,'Marker','o','MarkerFaceColor',c(i,:),'LineStyle','none','MarkerEdgeColor',c(i,:));
	L{i} = ['\tau_{gain} = ' oval(4*tau_gain(i)) 'ms'];
end 
lh4 = legend(l,L);
lh4.Position = [0.795   0.62    0.038    0.08];
xlabel(ax(7),'Gain control timescale (ms)')
ylabel(ax(7),['Correlation between' char(10) 'gain and \mu_{stimulus}'])

prettyFig('fs',14);

if being_published
	snapnow
	delete(gcf)
end

% move some plots
ax(3).Position(1) = .76;
ax(7).Position(1) = .76;



%% Version Info
%
pFooter;




