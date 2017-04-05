

pHeader;

% supp figure
supp_fig = figure('outerposition',[0 0 1500 750],'PaperUnits','points','PaperSize',[1500 750]); hold on
for i = 10:-1:1
	axs(i) = subplot(2,5,i); hold on
end

%  ######  #### ##     ## ##     ## ##          ###    ######## ####  #######  ##    ##  ######  
% ##    ##  ##  ###   ### ##     ## ##         ## ##      ##     ##  ##     ## ###   ## ##    ## 
% ##        ##  #### #### ##     ## ##        ##   ##     ##     ##  ##     ## ####  ## ##       
%  ######   ##  ## ### ## ##     ## ##       ##     ##    ##     ##  ##     ## ## ## ##  ######  
%       ##  ##  ##     ## ##     ## ##       #########    ##     ##  ##     ## ##  ####       ## 
% ##    ##  ##  ##     ## ##     ## ##       ##     ##    ##     ##  ##     ## ##   ### ##    ## 
%  ######  #### ##     ##  #######  ######## ##     ##    ##    ####  #######  ##    ##  ######  

% ##    ##    ###    ########     ######  ######## #### ##     ## 
% ###   ##   ## ##      ##       ##    ##    ##     ##  ###   ### 
% ####  ##  ##   ##     ##       ##          ##     ##  #### #### 
% ## ## ## ##     ##    ##        ######     ##     ##  ## ### ## 
% ##  #### #########    ##             ##    ##     ##  ##     ## 
% ##   ### ##     ##    ##       ##    ##    ##     ##  ##     ## 
% ##    ## ##     ##    ##        ######     ##    #### ##     ## 


cdata = consolidateData2('/local-data/DA-paper/data-for-paper/nat-stim/ab3-2ac');
v2struct(cdata)

% first, fit the DA model to the data. (this was pre-fit, we're simply loading the fit here)
load('/local-data/DA-paper/data-for-paper/fig7/DA_Model_fit_to_naturalistic_data.mat','p')


% generate DA model responses
dd = ORNData;
for i = 1:max(orn)
	S = PID(:,orn==i);
	dd(i).stimulus = nanmean(S,2);
	dd(i).stimulus = dd(i).stimulus - mean(dd(i).stimulus(1:5e3));
	dd(i).firing_rate = DAModelv2(dd(i).stimulus,p(i));
end
dd([1 4]) = [];
p([1 4]) = [];
dd = backOutFilters(dd);


history_lengths = round(logspace(1,log10(3e4),30)); % all the history lengths we look at, in ms
history_lengths(28) = [];
example_history_lengths = 300;
gain_mu = struct; gain_mu(length(example_history_lengths)).gain = []; gain_mu(1).mu = [];
rho = NaN(length(history_lengths),length(dd));

use_this_segment = false(length(dd(1).firing_rate),1);
use_this_segment(5e3:end-5e3) = true;


for i = 1:length(dd)
	temp = dd(i);
	pred = nanmean(temp.firing_projected,2); pred = pred(use_this_segment);
	resp = nanmean(temp.firing_rate,2);  resp = resp(use_this_segment);
	stim = nanmean(temp.stimulus,2); stim = stim - mean(stim(1:5e3));
	stim = stim(use_this_segment);

	% find when the valve opens
	[ons,offs] = findWhiffs(stim);
	ons = ons + finddelay(stim,resp);
	offs = offs + finddelay(stim,resp);

	% plot the gain in each of these windows
	[gain,gain_err] = findGainInWindows(ons,offs,pred,resp);

	rm_this = gain < 0 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];
	gain_err(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	for j = 1:length(example_history_lengths)
		mu = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_lengths(j)));
		gain_mu(j).mu = [gain_mu(j).mu(:); mu(:)];
		gain_mu(j).gain = [gain_mu(j).gain(:); gain(:)];
	end


	% also find rho for various values of the history length and plot it
	rho(:,i) = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

end

i = 1;
axes(axs(1))
plotPieceWiseLinear(gain_mu(i).mu,gain_mu(i).gain,'Color',[0 0 0]);
xlabel(axs(1),'\mu_{Stimulus} in preceding 300 ms (V)','Color','k')
ylabel(axs(1),'DA model gain (Hz/V)','Color','k')
set(axs(1),'XScale','log','YScale','log','XLim',[1e-4 11],'XTick',[1e-3 1e-2 1e-1 1e0 10])

% show rho vs. history length for the simulations
errorbar(axs(6),history_lengths,nanmean(rho,2),nanstd(rho,[],2)/sqrt(length(dd)),'k');
set(axs(6),'XScale','log','YLim',[-1 1],'XTick',[10 1e2 1e3 1e4],'XLim',[10 3e4])
xlabel(axs(6),'History length (ms)')
ylabel(axs(6),['Correlation between' char(10) 'gain and \mu_{stimulus}'])

% show the timescale of gain control of the model
for i = 1:length(p)
	plot(axs(6),[mean([p(i).tau_z].*[p(i).n_z]) mean([p(i).tau_z].*[p(i).n_z])],[-1 1],'r--')
end

title(axs(1),['Naturalistic stimulus' char(10) 'DA model simulations'])

clearvars -except axs being_published supp_fig lh

% ##     ##  ######   ######   
% ###   ### ##    ## ##    ##  
% #### #### ##       ##        
% ## ### ##  ######  ##   #### 
% ##     ##       ## ##    ##  
% ##     ## ##    ## ##    ##  
% ##     ##  ######   ######   

% define colors, etc. 
c = lines(10);
LFP_color = c(4,:);
firing_color = c(5,:);
model_color = [1 0 0];

% define what we want to work on
data_hashes = {'f11c4a5792d0c9fec7c40fd6aa2fce40','93ba5d68174e3df9f462a1fc48c581da','bcd4cf4fe12817d084a2b06f981161ee','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1'};
odour_names = {'isoamyl-acetate','ethyl-acetate','1-pentanol','1-pentanol','2-butanone'};
orn_names = {'pb1A','ab3A','ab3A','ab2A','ab2A',''};

markers = {'+','d','o','x','*'};
% core loop
for i = 1:length(data_hashes)
	clear cdata
	cdata = consolidateData2(getPath(dataManager,data_hashes{i}));
	if i < 4
		cdata.a = 25e3; cdata.z = 45e3;
	end
	if i < 2
		cdata.PID(:,mean(cdata.PID) > .35) = NaN;
	end


	cdata = cleanMSGdata(cdata,'extract_filter',false);

	plot_handles = plotMSGKinetics(cdata,axs(10));
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
for i = 1:length(markers)
	l(i) = plot(axs(10),NaN,NaN,'Marker',markers{i},'Color','k','LineStyle','none');
	L{i} = [orn_names{i} ' ' odour_names{i}];
end
lh(1) = legend(l,L,'Location','southeast');
lh(1).Position = [0.7860 0.6188 0.1528 0.2369];
lh(1).FontSize = 15;

set(axs(10),'YLim',[0 200],'XLim',[-0.1 1.1])
xlabel(axs(10),'Mean stimulus (rescaled)')


clearvars -except axs being_published supp_fig lh

axes(axs(5))
axis off
title(axs(5),['Gaussian stimulus' char(10) 'LFP and firing lags'])

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
	load(getPath(dataManager,'86946ed05ec73186d8371166583141ba'))
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

	rm_this =  gain < 0 | gain_err < .8;
	gain(rm_this) = [];
	ons(rm_this) = [];
	offs(rm_this) = [];

	% find the mean stimulus in the preceding X ms in these windows
	mean_stim = findMeanInWindows(ons,offs,computeSmoothedStimulus(stim,example_history_length));

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);

	plot(axs(9),history_lengths,rho,'-','Color',c(i,:));
end

set(axs(9),'XScale','log','YLim',[-1 0],'XTick',[1e1 1e2 1e3 1e4],'XLim',[10 1e4],'XLim',[10 1.1e4])
xlabel(axs(9),'History length (ms)')

% fake some plots for a nice legend
clear L l
for i = 1:length(do_these)
	l(i) = plot(axs(6),NaN,NaN,'Marker','o','MarkerFaceColor',c(i,:),'LineStyle','none');
	L{i} = [orn_data(do_these(i)).neuron_name ' ' odour_names{i}];
end
lh(2) = legend(l,L,'Location','southeast');
lh(2).Position = [0.6200 0.6367 0.1500 0.2000];
lh(2).FontSize = 15;

axes(axs(4))
axis off
title(axs(4),['Binary stimulus' char(10) 'Gain control timescale'])

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


%  ######  #### ##     ## ##     ## ##          ###    ######## ####  #######  ##    ##  ######  
% ##    ##  ##  ###   ### ##     ## ##         ## ##      ##     ##  ##     ## ###   ## ##    ## 
% ##        ##  #### #### ##     ## ##        ##   ##     ##     ##  ##     ## ####  ## ##       
%  ######   ##  ## ### ## ##     ## ##       ##     ##    ##     ##  ##     ## ## ## ##  ######  
%       ##  ##  ##     ## ##     ## ##       #########    ##     ##  ##     ## ##  ####       ## 
% ##    ##  ##  ##     ## ##     ## ##       ##     ##    ##     ##  ##     ## ##   ### ##    ## 
%  ######  #### ##     ##  #######  ######## ##     ##    ##    ####  #######  ##    ##  ######  

% first, fit a DA model to the methyl butyrate ab3A data
clear p
p.   s0 = -8.5586e-04;
p.  n_z = 2;
p.tau_z = 48.3750;
p.  n_y = 1.9688;
p.tau_y = 15.6406;
p.    C = 0.0126;
p.    A = 2.5139e+04;
p.    B = 423.0312;

% generate synthetic data responses and do the analyses on that

% tweak the model's tau_gain
tau_gain = [50 100 200 400];
actual_tau_gain = NaN*tau_gain;

resp = mean(orn_data(16).firing_rate,2);
resp = resp(30e3:end);


% return

c = parula(length(tau_gain)+1);
clear l
for i = 1:length(tau_gain)
	textbar(i,length(tau_gain));

	stim = nanmean(orn_data(16).stimulus,2);
	% stim = stim/nanmean(stim);
	% generate responses
	p.tau_z = tau_gain(i);
	[resp,~,~,~,Kz] = DAModelv2(stim,p);
	[~,actual_tau_gain(i)] = max(Kz);

	resp = resp(30e3:end);
	stim = stim(30e3:end);


	% fit a filter to this
	K = fitFilter2Data(stim,resp,'filter_length',1e3,'offset',200);
	K = K(100:end-100);
	filtertime = (1:length(K)) - 100;
	pred = convolve(1:length(stim),stim,K,filtertime);


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
	l(i) = plot(axs(2),mean_stim,gain(idx),'+-','Color',c(i,:));

	% also find rho for various values of the history length and plot it
	rho = findRhoForHistoryLengths(gain,stim,ons,offs,history_lengths);
	[~,temp] = min(rho);
	min_hist_lengths(i) = history_lengths(temp);


	plot(axs(7),history_lengths,rho,'.-','Color',c(i,:))
end
set(axs(7),'XScale','log','YLim',[-1 0],'XTick',[10 100 1e3 1e4],'XLim',[10 1.1e4])
set(axs(2),'XScale','log','YScale','log','XLim',[5e-4 5e-3],'YLim',[.1 10],'XTick',[1e-3 2e-3 4e-3])
title(axs(2),['Binary stimulus' char(10) 'DA model simulations'])
xlabel(axs(2),'\mu_{Stimulus} in preceding 300 ms (V)')

clear l L
for i = 1:length(tau_gain)
	l(i) = plot(axs(2),NaN,NaN,'Marker','o','MarkerFaceColor',c(i,:),'LineStyle','none','MarkerEdgeColor',c(i,:));
	L{i} = ['\tau_{gain} = ' oval(p.n_z*tau_gain(i)) 'ms'];
end 
lh(4) = legend(l,L);
lh(4).Position = [0.3121 0.6241 0.0736 0.0952];
xlabel(axs(7),'History length (ms)')


prettyFig(supp_fig,'fs',14);

% make some room between the top and bottom row
for i = 1:5
	axs(i).Position(4) = .33;
	axs(i).Position(2) = .6;
end

for i = 6:10
	axs(i).Position(4) = .33;
end

% fix some positions

axs(1).Position(1) = axs(1).Position(1) - .01;
axs(2).Position(1) = axs(2).Position(1) - .01;

axs(6).Position(1) = axs(1).Position(1);
axs(6).Position(3) = axs(1).Position(3);

axs(7).Position(1) = axs(2).Position(1);
axs(7).Position(3) = axs(2).Position(3);


% deintersect some axes
for i = 1:10
	try
		deintersectAxes(axs(i))
	catch
	end
end

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;
