% fig_kinetics.m
% 
% created by Srinivas Gorur-Shandilya at 2:14 , 07 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;
dm = dataManager;

figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on
for i = 6:-1:1
	ax(i) = subplot(2,3,i); hold on
end
delete(ax(3)); delete(ax(6))


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

	axes(ax(4))

	plotPieceWiseLinear(t,lag,'Color',c(1,:),'nbins',19);


	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,X,1e3,25);
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
	lag(rm_this) = [];
	mean_x(rm_this) = [];
	t = time_since_thresh_crossing;
	t(rm_this) = []; t(t<10) = NaN;

	axes(ax(1))
	l(1) = plotPieceWiseLinear(mean_x,lag,'Color',c(2,:),'nbins',19);

	axes(ax(4))
	plotPieceWiseLinear(t,lag,'Color',c(2,:),'nbins',19);

end
xlabel(ax(1),'\mu_{Stimulus} in preceding 1s (V)')
ylabel(ax(1),'Lag (ms)')
set(ax(1),'YLim',[0 140],'XLim',[0 0.6])
L = legend(l,{'LFP','Firing Rate'},'Location','southeast');
title(ax(1),['ab3A' char(10) 'ethyl-acetate'])

set(ax(4),'YLim',[0 140],'XLim',[10 5000],'XScale','log')
xlabel(ax(4),'Time since odor encounter (ms)')
ylabel(ax(4),'Lag (ms)')


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

	plot_handles = plotMSGKinetics(cdata,ax(2));

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
	%title(t);
end

% fake some plots for a nice legend
clear l L
for i = 1:length(markers)-1
	l(i) = plot(ax(2),NaN,NaN,'Marker',markers{i},'Color','k','LineStyle','none');
	L{i} = [orn_names{i} ' ' odour_names{i}];
end
lh = legend(l,L,'Location','southeast');
lh.Position = [0.66 0.7 0.15 0.1];
lh.FontSize = 15;

set(ax(2),'YLim',[0 200],'XLim',[-0.1 1.1])
xlabel(ax(2),'Mean stimulus (rescaled)')


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

clear plot_handles
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

	plot_handles(i) = plot(ax(5),history_lengths,rho,'-');
end

set(ax(5),'XScale','log','YLim',[-1 0],'XTick',[1e2 1e3 1e4],'XLim',[100 1e4])
xlabel(ax(5),'Gain control timescale (ms)')
ylabel(ax(5),'Spearman''s \rho')


% fake some plots for a nice legend
clear L
for i = 1:length(plot_handles)
	L{i} = [orn_data(do_these(i)).neuron_name ' ' odour_names{i}];
end
lh = legend(plot_handles,L,'Location','southeast');
lh.Position = [0.66 0.15 0.15 0.2];
lh.FontSize = 15;

prettyFig('fs',.5,'font_units','centimeters')

% deintersect the bottom row
for i = 4:5
	deintersectAxes(ax(i));
end

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


