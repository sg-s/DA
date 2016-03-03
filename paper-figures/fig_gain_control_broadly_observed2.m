% fast gain control is broadly observed 
% makes figure 6 of the paper, showing that gain control is broadly observed
% 
% created by Srinivas Gorur-Shandilya at 9:56 , 04 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

%      ########     ###             ##      ## #### ########  ######## ##       ##    ## 
%      ##     ##   ## ##            ##  ##  ##  ##  ##     ## ##       ##        ##  ##  
%      ##     ##  ##   ##           ##  ##  ##  ##  ##     ## ##       ##         ####   
%      ##     ## ##     ##          ##  ##  ##  ##  ##     ## ######   ##          ##    
%      ##     ## #########          ##  ##  ##  ##  ##     ## ##       ##          ##    
%      ##     ## ##     ##          ##  ##  ##  ##  ##     ## ##       ##          ##    
%      ########  ##     ##           ###  ###  #### ########  ######## ########    ## 

%       #######  ########   ######  ######## ########  ##     ## ######## ########  
%      ##     ## ##     ## ##    ## ##       ##     ## ##     ## ##       ##     ## 
%      ##     ## ##     ## ##       ##       ##     ## ##     ## ##       ##     ## 
%      ##     ## ########   ######  ######   ########  ##     ## ######   ##     ## 
%      ##     ## ##     ##       ## ##       ##   ##    ##   ##  ##       ##     ## 
%      ##     ## ##     ## ##    ## ##       ##    ##    ## ##   ##       ##     ## 
%       #######  ########   ######  ######## ##     ##    ###    ######## ########  

% load the data
if ~exist('orn_data','var')
	load('Carlotta_Data.mat')
end


%% Detailed Report: Experimental Replicates
% In this section, we look at all the data from Carlotta in detail, ORN-by-ORN. First, we start off by looking at her experimental replicates, where the ab3A ORN responds to methyl butyrate. 

% first column: experimental replicates
do_these = [7 9 13 15 16];

for i = 1:length(do_these)
	figure('outerposition',[0 0 1300 800],'PaperUnits','points','PaperSize',[1300 800]); hold on
	for j = 1:5
		ax(j) = subplot(2,3,j); hold on
	end

	% show the filter
	plot(orn_data(do_these(i)),ax(1),'Filter.firing_rate');

	% show the nonlinearity -- raw 
	plot(orn_data(do_these(i)),ax(2),'ioCurve.firing_rate','show_NL',false,'nbins',50);

	% now fit a NL
	temp = orn_data(do_these(i));
	temp = fitPieceWiseNL(temp);
	plot(temp,ax(3:5),'valveGainAnalysis.firing_rate.mu','history_lengths',logspace(1,4,30));

	% cosmetics
	set(ax(4),'YLim',[.5 2],'YScale','log')

	prettyFig('fs=18;')
	if being_published
		snapnow
		delete(gcf)
	end

end

%% Detailed Report: Different Odours
% We now see if this is true for different odours. 

do_these = [7 8 10 14 17 18];

for i = 1:length(do_these)
	figure('outerposition',[0 0 1300 800],'PaperUnits','points','PaperSize',[1300 800]); hold on
	for j = 1:5
		ax(j) = subplot(2,3,j); hold on
	end

	% show the filter
	plot(orn_data(do_these(i)),ax(1),'Filter.firing_rate');

	% show the nonlinearity -- raw 
	plot(orn_data(do_these(i)),ax(2),'ioCurve.firing_rate','show_NL',false,'nbins',50);

	% now fit a NL
	temp = orn_data(do_these(i));
	temp = fitPieceWiseNL(temp);
	plot(temp,ax(3:5),'valveGainAnalysis.firing_rate.mu','history_lengths',logspace(1,4,30));

	% cosmetics
	set(ax(4),'YLim',[.5 2],'YScale','log')

	prettyFig('fs=18;')

	suptitle(['Odour: ' orn_data(do_these(i)).odour_name])
	if being_published
		snapnow
		delete(gcf)
	end

end


%% Detailed Analysis: Different ORNs
% We now compare fast gain control in different ORNs

do_these = [11 17 25];



% make the figures
figure('outerposition',[0 0 1300 800],'PaperUnits','points','PaperSize',[1300 800]); hold on
clear axes_handles
for i = 1:12
	axes_handles(i) = subplot(3,4,i);
	hold on
end


title(axes_handles(1),'Experimental Replicates')
title(axes_handles(2),'Different Odors')
title(axes_handles(3),'Different ORNs')
title(axes_handles(4),'Locust EAG')


clear plot_handles
for i = 1:length(do_these)
	do_this = do_these(i);
	plot_handles(i).h = plot(orn_data(do_this),axes_handles([5 9]),'instGainAnalysis.firing_rate.mu','history_lengths',round(logspace(2,4,50)),'data_bin_type','dots','history_length',200);
end


% cosmetics
c =  lines(length(do_these));
for i = 1:length(do_these)
	set(plot_handles(i).h(1).lines,'Color',c(i,:),'LineStyle','none','Marker','+')
	set(plot_handles(i).h(2).lines,'Color',c(i,:),'LineStyle','none','Marker','+')
end

set(axes_handles(5),'YLim',[0.5 2],'XLim',[0 0.9])
set(axes_handles(9),'YLim',[-1 0])

% show a explanatory graphic
o = imread('../images/exp-rep.png');
axes(axes_handles(1));
imagesc(o);
axis ij
axis image
axis off
uistack(axes_handles(1),'bottom')


% get the legend right
ph = [plot_handles.h]; ph = ph(1:2:end);
L = {'06/03','05/28','06/05','06/12','06/19'};
lh = legend(ph,L);
set(lh,'Position',[0.2292 0.75 0.0500 0.0770])

%     ########  #### ######## ######## ######## ########  ######## ##    ## ######## 
%     ##     ##  ##  ##       ##       ##       ##     ## ##       ###   ##    ##    
%     ##     ##  ##  ##       ##       ##       ##     ## ##       ####  ##    ##    
%     ##     ##  ##  ######   ######   ######   ########  ######   ## ## ##    ##    
%     ##     ##  ##  ##       ##       ##       ##   ##   ##       ##  ####    ##    
%     ##     ##  ##  ##       ##       ##       ##    ##  ##       ##   ###    ##    
%     ########  #### ##       ##       ######## ##     ## ######## ##    ##    ##    
    
%      #######  ########   #######  ##     ## ########   ######  
%     ##     ## ##     ## ##     ## ##     ## ##     ## ##    ## 
%     ##     ## ##     ## ##     ## ##     ## ##     ## ##       
%     ##     ## ##     ## ##     ## ##     ## ########   ######  
%     ##     ## ##     ## ##     ## ##     ## ##   ##         ## 
%     ##     ## ##     ## ##     ## ##     ## ##    ##  ##    ## 
%      #######  ########   #######   #######  ##     ##  ######  


do_these = [7 8 10 14 17 18];

clear plot_handles
for i = 1:length(do_these)
	do_this = do_these(i);
	plot_handles(i).h = plot(orn_data(do_this),axes_handles([6 10]),'gain_analysis.binned','mean_stim_bins',200);
end

% replot the binned gain vs. mean stim with futher binning
for i = 1:length(do_these)
	x = get(plot_handles(i).h(1),'XData');
	y = get(plot_handles(i).h(1),'YData');
	rm_this = x < .1 | x > .9;
	x(rm_this) = []; y(rm_this) = [];
	delete(plot_handles(i).h(1));
	axes(axes_handles(6));
	plot_handles(i).h(1) = plotPieceWiseLinear(x,y,'nbins',5,'use_std',true);
end

% cosmetics
c =  lines(length(do_these));
for i = 1:length(do_these)
	set(plot_handles(i).h(1),'Color',c(i,:),'LineStyle','-','Marker','none')
	set(plot_handles(i).h(2),'Color',c(i,:),'LineStyle','-','Marker','none')
end

set(axes_handles(6),'YLim',[0.5 2],'XLim',[0 0.9])
set(axes_handles(10),'YLim',[-1 0])

% show a explanatory graphic
o = imread('../images/diff-odors.png');
axes(axes_handles(2));
imagesc(o);
axis ij
axis image
axis off
uistack(axes_handles(2),'bottom')


% get the legend right
ph = [plot_handles.h]; ph = ph(1:2:end);
odours = {'1but','1o3ol','dsucc','2ac','2but','5ol'};
lh = legend(ph,odours);
set(lh,'Position',[0.35 0.75 0.0300 0.0770],'Box','off')


%    ########  #### ######## ######## ######## ########  ######## ##    ## ######## 
%    ##     ##  ##  ##       ##       ##       ##     ## ##       ###   ##    ##    
%    ##     ##  ##  ##       ##       ##       ##     ## ##       ####  ##    ##    
%    ##     ##  ##  ######   ######   ######   ########  ######   ## ## ##    ##    
%    ##     ##  ##  ##       ##       ##       ##   ##   ##       ##  ####    ##    
%    ##     ##  ##  ##       ##       ##       ##    ##  ##       ##   ###    ##    
%    ########  #### ##       ##       ######## ##     ## ######## ##    ##    ##    
   
%    ##    ## ######## ##     ## ########   #######  ##    ##  ######  
%    ###   ## ##       ##     ## ##     ## ##     ## ###   ## ##    ## 
%    ####  ## ##       ##     ## ##     ## ##     ## ####  ## ##       
%    ## ## ## ######   ##     ## ########  ##     ## ## ## ##  ######  
%    ##  #### ##       ##     ## ##   ##   ##     ## ##  ####       ## 
%    ##   ### ##       ##     ## ##    ##  ##     ## ##   ### ##    ## 
%    ##    ## ########  #######  ##     ##  #######  ##    ##  ######  

do_these = [11 17 25];

clear plot_handles
for i = 1:length(do_these)
	do_this = do_these(i);
	plot_handles(i).h = plot(orn_data(do_this),axes_handles([7 11]),'gain_analysis.binned','mean_stim_bins',200);
end

% replot the binned gain vs. mean stim with futher binning
for i = 1:length(do_these)
	x = get(plot_handles(i).h(1),'XData');
	y = get(plot_handles(i).h(1),'YData');
	rm_this = x < .1 | x > .9;
	x(rm_this) = []; y(rm_this) = [];
	delete(plot_handles(i).h(1));
	axes(axes_handles(7));
	plot_handles(i).h(1) = plotPieceWiseLinear(x,y,'nbins',5,'use_std',true);
end

% cosmetics
c =  lines(length(do_these));
for i = 1:length(do_these)
	set(plot_handles(i).h(1),'Color',c(i,:),'LineStyle','-','Marker','none')
	set(plot_handles(i).h(2),'Color',c(i,:),'LineStyle','-','Marker','none')
end

set(axes_handles(7),'YLim',[0.5 2],'XLim',[0 0.9])
set(axes_handles(11),'YLim',[-1 0])

% show a explanatory graphic
o = imread('../images/diff-orns.png');
axes(axes_handles(3));
imagesc(o);
axis ij
axis image
axis off
uistack(axes_handles(3),'bottom')


% get the legend right
ph = [plot_handles.h]; ph = ph(1:2:end);
lh = legend(ph,odours,{'','',''});
set(lh,'Position',[15 0.75 0.0300 0.0770],'Box','off')


% ##        #######   ######  ##     ##  ######  ######## 
% ##       ##     ## ##    ## ##     ## ##    ##    ##    
% ##       ##     ## ##       ##     ## ##          ##    
% ##       ##     ## ##       ##     ##  ######     ##    
% ##       ##     ## ##       ##     ##       ##    ##    
% ##       ##     ## ##    ## ##     ## ##    ##    ##    
% ########  #######   ######   #######   ######     ##    

% load data
if ~exist('locust_data','var')
	load('/local-data/DA-paper/locust/example-data')

	% clean up, sub-sample to 1ms
	PID = PID1; clear PID1
	EAG = EAG1; clear EAG1 

	PID = PID(:,1:10:end)';
	EAG = EAG(:,1:10:end)';
	valve = ODR1(:,1:10:end)';
	valve(valve<max(max(valve))/2) = 0;
	valve(valve>0) = 1;

	% set zero
	for i = 1:width(PID)
		PID(:,i) = PID(:,i) - mean(PID(1:300,i));
		EAG(:,i) = EAG(:,i) - mean(EAG(1:300,i));
		% filter
		PID(:,i) = bandPass(PID(:,i),Inf,30);
		EAG(:,i) = bandPass(EAG(:,i),2e3,Inf);
	end


	locust_data = ORNData;
	locust_data.filtertime_LFP = -.5:1e-3:.9;
	locust_data.regularisation_factor = 1;
	locust_data.stimulus = PID;
	locust_data.LFP = EAG;
	locust_data = backOutFilters(locust_data);
end

clear plot_handles
plot_handles.h = plot(locust_data,axes_handles([8 12]),'gain_analysis.LFP','mean_stim_bins',200,'history_length',1);


% replot the binned gain vs. mean stim with futher binning
x = get(plot_handles.h(1),'XData');
y = get(plot_handles.h(1),'YData');
rm_this = y < 0;
x(rm_this) = []; y(rm_this) = [];
delete(plot_handles.h(1));
axes(axes_handles(8));
plot_handles.h(1) = plotPieceWiseLinear(x,y,'nbins',5,'use_std',true);


set(plot_handles.h(1),'LineStyle','-','Marker','none')
set(plot_handles.h(2),'LineStyle','-','Marker','none')

set(axes_handles(12),'YLim',[-1 0])

prettyFig('plw=1.5;','lw=1.5;','fs=14;')

% turn all boxes off
for i = 1:length(axes_handles)
	set(axes_handles(i),'box','off')
end


if being_published
	snapnow
	delete(gcf)
end


%% Version Info
% The file that generated this document is called:
pFooter;
