% fig_kinetics.m
% 
% created by Srinivas Gorur-Shandilya at 2:14 , 07 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;
dm = dataManager;

figure('outerposition',[0 0 1400 600],'PaperUnits','points','PaperSize',[1400 600]); hold on

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

top_row(1) = subplot(2,6,1); hold on

c = lines(3);

min_acceptable_corr = .5;
min_acceptable_lag = 2;
clear l
for i = 1:length(od)
	S = nanmean(od(i).stimulus,2); 
	R = nanmean(od(i).firing_rate,2);
	X = -nanmean(od(i).LFP,2);

	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,R,1e3,25);
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
	lag(rm_this) = [];
	mean_x(rm_this) = [];

	l(2) = plotPieceWiseLinear(mean_x,lag,'Color',c(1,:),'nbins',19);


	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,X,1e3,25);
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
	lag(rm_this) = [];
	mean_x(rm_this) = [];

	l(1) = plotPieceWiseLinear(mean_x,lag,'Color',c(2,:),'nbins',19);
end
xlabel('\mu_{Stimulus} in preceding 1s (V)')
ylabel('Lag (ms)')
set(gca,'YLim',[0 200],'XLim',[0 0.6])
title(['ab3A' char(10) 'ethyl-acetate'])
L = legend(l,'LFP','Firing Rate');

% ##     ##  ######   ######   
% ###   ### ##    ## ##    ##  
% #### #### ##       ##        
% ## ### ##  ######  ##   #### 
% ##     ##       ## ##    ##  
% ##     ## ##    ## ##    ##  
% ##     ##  ######   ######   


% define what we want to work on
data_hashes = {'28ee201995ff7a193f072bd30f556348','d06180a4de81dea98526d349e925dd40','dbaf184b0b69939fb4d949eb4bd38995','ae5f11c1671c863a09c4bf9ec683ec16','a9ce591b1ec29a56c35db8f19b374b97'};
odour_names = {'ethyl-acetate','1-pentanol','1-pentanol','2-butanone','isoamyl-acetate'};
orn_names = {'ab3A','ab3A','ab2A','ab2A','pb1A'};


% core loop
for i = 1:length(data_hashes)
	clear cdata

	cdata = consolidateData2(dm.getPath(data_hashes{i}));
	if i < 3
		cdata.a = 25e3; cdata.z = 45e3;
	end
	cdata = cleanMSGdata(cdata,'extract_filter',false);

	top_row(i+1) = subplot(2,length(data_hashes)+1,i+1); hold on
	plotMSGKinetics(cdata,top_row(i+1));
	set(gca,'YLim',[0 200],'XLim',[0 max(nanmean(cdata.PID))*1.1])
	t = [orn_names{i} char(10) odour_names{i}];
	title(t);

end

%% Now, show the timescale of gain control using Carlotta's data

%% global parameters
history_lengths = round(logspace(1.7,4,50)); % all the history lengths we look at, in ms
example_history_length = 300; % this history length shown in the first row, in ms

% load the data
if ~exist('orn_data','var')
	load('/local-data/DA-paper/fig4/Carlotta_Data.mat')
end

do_these = [18 7 8 10 14 17 12];
odour_names = {'1-pentanol','methyl-butyrate','1-octen-3-ol','diethyl-succinate','ethyl-acetate','2-butanone','isoamyl-acetate'};

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

	bottom_row(i) = subplot(2,7,7+i); hold on
	plot(bottom_row(i),history_lengths,rho,'k.-')
	set(bottom_row(i),'XScale','log','YLim',[-1 0],'XTick',[1e2 1e3 1e4],'XLim',[100 1e4])
	xlabel(bottom_row(i),'Gain control timescale (ms)')
	ylabel(bottom_row(i),'Spearman''s \rho')
	t = [orn_data(do_these(i)).neuron_name char(10) odour_names{i}];
	title(bottom_row(i),t);
	drawnow;

end

% cosmetic fixes:

% make top row smaller, move them up a bit
for i = 1:length(top_row)
	top_row(i).Position(2) = .6;
	top_row(i).Position(4) = .3;
end

% remove extra labels from top row
for i = 2:length(top_row)
	ylabel(top_row(i),'');
end

% spread out the bottom row a bit
padding_left = .07;
padding_right = 0;
for i = 1:length(bottom_row)
	bottom_row(i).Position(1) = (1-padding_left-padding_right)*(i-1)*(1/length(bottom_row)) + padding_left;
	bottom_row(i).Position(4) = .3;
end

prettyFig('fs',.5,'font_units','centimeters')

% deintersect the bottom row
for i = 1:length(bottom_row)
	deintersectAxes(bottom_row(i));
end

L.Box = 'off';

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


