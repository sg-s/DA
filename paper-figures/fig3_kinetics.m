% fig_kinetics.m
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
LFP_firing_lag_plot = ax(4);
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
temp = od([1 3]); % other entries are blank, because we have no data. 
clear od;
load('/local-data/DA-paper/data-for-paper/nat-stim/ab3A_nat_stim.ORNData','-mat')
od = [temp od];
clear temp


l_s = {'-','-','--','--','--','--','--','--','--','--'};
min_acceptable_corr = .5;
min_acceptable_lag = 5;
max_acceptable_lag = 300;
window_size = 1e3;
nbins = 40;
history_length = 300;
stim_thresh = .035;
clear l


for i = 1:length(od)
	textbar(i,length(od))
	S = nanmean(od(i).stimulus,2); S = S - mean(S(1:5e3));
	R = nanmean(od(i).firing_rate,2);
	X = -nanmean(od(i).LFP,2);

	% stim --> firing rate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	[lag, ~, max_corr] = findLagAndMeanInWindow(S,R,window_size,5);
	mean_x = vectorise(computeSmoothedStimulus(S,window_size));


	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];

	% then throw out some shitty data
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr | lag > max_acceptable_lag;
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];

	% plot
	[~,data] = plotPieceWiseLinear(log(mean_x),lag,'nbins',10,'make_plot',false,'proportional_bins',false);
	plot(lag_vs_stim_plot,exp(data.x),data.y,'Color',firing_color,'LineStyle',l_s{i});


 	% stm --> LFP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	[lag, ~, max_corr] = findLagAndMeanInWindow(S,X,window_size,5);
	mean_x = vectorise(computeSmoothedStimulus(S,window_size));

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];


	% then throw out some shitty data
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr | lag > max_acceptable_lag | mean_x < 1e-3;

	[~,data] = plotPieceWiseLinear(log(mean_x(~rm_this)),lag(~rm_this),'nbins',10,'make_plot',false,'proportional_bins',false);
	plot(lag_vs_stim_plot,exp(data.x),data.y,'LineStyle',l_s{i},'Color',LFP_color)


	% LFP --> fring rate

	[lag, ~, max_corr] = findLagAndMeanInWindow(X,R,window_size,5);
	mean_x = vectorise(computeSmoothedStimulus(S,window_size));
	mean_x = mean_x - min(mean_x);

	% first strip out the NaNs
	rm_this = isnan(lag);
	lag(rm_this) = []; mean_x(rm_this) = []; max_corr(rm_this) = [];


	% then throw out some shitty data
	rm_this = max_corr < min_acceptable_corr | lag < -max_acceptable_lag | lag > max_acceptable_lag | mean_x < 1e-3;

	[~,data] = plotPieceWiseLinear(log(mean_x(~rm_this)),lag(~rm_this),'nbins',10,'make_plot',false,'proportional_bins',false);
	plot(LFP_firing_lag_plot,exp(data.x),data.y,'LineStyle',l_s{i},'Color','k')


end

% labels -- main figure
xlabel(lag_vs_stim_plot,'\mu_{Stimulus} in preceding 300 ms (V)')
ylabel(lag_vs_stim_plot,'Lag (ms)')
set(lag_vs_stim_plot,'YLim',[0 160],'XLim',[1e-3 3],'XScale','log','YScale','linear','XTick',[1e-4 1e-3 1e-2 1e-1 1],'XMinorTick','off')



set(LFP_firing_lag_plot,'YLim',[-100 60],'XLim',[1e-3 3],'XScale','log','YScale','linear','XTick',[1e-4 1e-3 1e-2 1e-1 1],'XMinorTick','off')
xlabel(LFP_firing_lag_plot,'\mu_{Stimulus} in preceding 300 ms (V)')
ylabel(LFP_firing_lag_plot,'LFP \rightarrow firing lag (ms)')


% fake a fake legend for visual niceness
clear l
l(1) = plot(lag_vs_stim_plot,NaN,NaN,'Color',LFP_color);
l(2) = plot(lag_vs_stim_plot,NaN,NaN,'Color',firing_color);
L1 = legend(l,{'LFP','Firing rate'},'Location','southeast');

clear l
l(1) = plot(LFP_firing_lag_plot,NaN,NaN,'k-');
l(2) = plot(LFP_firing_lag_plot,NaN,NaN,'k--');
L2 = legend(l,{'ab2A','ab3A'},'Location','northwest');
L2.Position = [0.1700 0.4023 0.0800 0.0337];
uistack(LFP_firing_lag_plot,'top')


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
xlabel(rho_vs_history_length_plot,'History length (ms)')
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


