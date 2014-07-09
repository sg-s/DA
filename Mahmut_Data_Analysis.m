% Mahmut_Data_Analyis.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% specify which data to look at
clearvars -except options
load('/local-data/DA-paper/data.mat')
td = 7;
redo_bootstrap = 0;
whiff_anal = 0;

%% 
% How do ORNs respond to non-Gaussian inputs? Real odor stimuli are characterised by long tails and non-gaussian statisitcs, with large whiffs of odor that occur in periods of relatively low signal. Such a stimulus has been generated here, and the responses of ORNs to these signals is analysed in this figure.





%                     ########     ###    ########    ###    
%                     ##     ##   ## ##      ##      ## ##   
%                     ##     ##  ##   ##     ##     ##   ##  
%                     ##     ## ##     ##    ##    ##     ## 
%                     ##     ## #########    ##    ######### 
%                     ##     ## ##     ##    ##    ##     ## 
%                     ########  ##     ##    ##    ##     ## 

%% Rough Overview of Data
% The following figure shows what the stimulus and the neuron's response looks like. 

% first, remove the baseline from the PID
if isvector(data(td).PID)
	data(td).PID = data(td).PID - mean(data(td).PID(1:3000));
	data(td).PID(data(td).PID < 0) = 0;
else
	% we have many PID trials
	for i = 1:width(data(td).PID)
		data(td).PID(:,i) = data(td).PID(:,i) - mean(data(td).PID(1:3000,i));
		data(td).PID(data(td).PID(:,i) < 0,i) = 0;
	end
end


figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(4,1,1), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
titlestr = strkat(data(td).neuron,' response to ',data(td).odor);
title(titlestr)
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
set(gca,'YLim',[-0.1 1.1*max(max(data(td).PID))])

subplot(4,1,2), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])
a = min(nonzeros(data(td).PID(:)));
z = 1.1*max((data(td).PID(:)));
set(gca,'YTick',logspace(log10(a),log10(z),4));
set(gca,'YScale','log','YMinorTick','on')

subplot(4,1,3), hold on
if isfield(data,'spiketimes') && isfield(data,'spiketimeB')
	if isempty(data(td).spiketimeB)
		raster2(data(td).spiketimes)
	else
		raster2(data(td).spiketimeB,data(td).spiketimes)
	end
end
set(gca,'XLim',[min(data(td).time) max(data(td).time)])

subplot(4,1,4), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
set(gca,'YLim',[0 1.1*max(max(data(td).ORN))])

PrettyFig;

snapnow;
delete(gcf);


%     ##     ##    ###    ########  ####    ###    ########  ##       #### ######## ##    ## 
%     ##     ##   ## ##   ##     ##  ##    ## ##   ##     ## ##        ##     ##     ##  ##  
%     ##     ##  ##   ##  ##     ##  ##   ##   ##  ##     ## ##        ##     ##      ####   
%     ##     ## ##     ## ########   ##  ##     ## ########  ##        ##     ##       ##    
%      ##   ##  ######### ##   ##    ##  ######### ##     ## ##        ##     ##       ##    
%       ## ##   ##     ## ##    ##   ##  ##     ## ##     ## ##        ##     ##       ##    
%        ###    ##     ## ##     ## #### ##     ## ########  ######## ####    ##       ##    

%% Trial to Trial variability
% How variable is the response and the stimulus from trial to trial? The following figure shows individual traces of the stimulus and the neuron's response for each trial. 


figure('outerposition',[0 0 1000 1000],'PaperUnits','points','PaperSize',[1000 1000]); hold on
set(gcf,'DefaultAxesColorOrder',jet(width(data(td).PID)))
subplot(3,2,1), hold on
plot(data(td).time,data(td).PID);
ylabel('PID (V)')
set(gca,'XLim',[25 35])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

subplot(3,2,2), hold on
plot(data(td).time,data(td).PID);
ylabel('PID (V)')
set(gca,'XLim',[38 48])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

if isfield(data,'spiketimes') && isfield(data,'spiketimeB')
	subplot(3,2,3), hold on
	raster2(data(td).spiketimes)
	set(gca,'XLim',[25 35])


	subplot(3,2,4), hold on
	raster2(data(td).spiketimes)
	set(gca,'XLim',[38 48])
end

subplot(3,2,5), hold on
plot(data(td).time,data(td).ORN);
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[25 35])
set(gca,'YLim',[0 1.1*max(max(data(td).ORN))])

subplot(3,2,6), hold on
plot(data(td).time,data(td).ORN);
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[38 48])
set(gca,'YLim',[0 1.1*max(max(data(td).ORN))])
PrettyFig;

snapnow;
delete(gcf);

%    ##      ## ##     ## #### ######## ########          ###    ##    ##    ###    ##       
%    ##  ##  ## ##     ##  ##  ##       ##               ## ##   ###   ##   ## ##   ##       
%    ##  ##  ## ##     ##  ##  ##       ##              ##   ##  ####  ##  ##   ##  ##       
%    ##  ##  ## #########  ##  ######   ######         ##     ## ## ## ## ##     ## ##       
%    ##  ##  ## ##     ##  ##  ##       ##             ######### ##  #### ######### ##       
%    ##  ##  ## ##     ##  ##  ##       ##             ##     ## ##   ### ##     ## ##       
%     ###  ###  ##     ## #### ##       ##             ##     ## ##    ## ##     ## ######## 

%%
% It looks like the amplitude of the signal is dropping every trial, but that the response amplitude isn't changing much (see the whiff at $t=39s$). 

%%
% To quantify this, let's attempt to find all the times when the stimulus is high, i.e., ten standard deviations above the baseline. This threshold ensures that we pick up all the whiffs, but ignore the blanks in between, as shown in the figure below. 

if whiff_anal

	figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
	plot(data(td).time, data(td).PID(:,1),'k'), hold on

	hs = data(td).PID(:,1) > 10*std(data(td).PID(1:3000,1));
	plot(data(td).time(hs), data(td).PID(hs,1),'.r','MarkerSize',24), hold on
	set(gca,'XLim',[28 36])
	xlabel('Time (s)')
	ylabel('Stimulus (V)')
	title('Whiff identification')
	set(gca,'YScale','log')
	PrettyFig;
	snapnow;
	delete(gcf);

end



%%
% Now we break up the trace so that we can perform a whiff-by-whiff analysis of the stimulus and the response, for each trial. The following figure shows how the whiff amplitude and response amplitude drop as a function of trial number.


if whiff_anal

	[whiff_ons,whiff_offs] = ComputeOnsOffs(hs);

	whiff_durations = whiff_offs - whiff_ons;
	whiff_ons(whiff_durations <20) = [];
	whiff_offs(whiff_durations <20) = [];
	whiff_durations = whiff_offs - whiff_ons;

	whiff_stim_max  = NaN(width(data(td).PID),length(whiff_ons));
	whiff_resp_max  = NaN(width(data(td).PID),length(whiff_ons));


	for j = 1:width(data(td).PID)
		for i = 1:length(whiff_ons)
			whiff_stim_max(j,i) = max(data(td).PID(whiff_ons(i):whiff_offs(i),j));
			whiff_resp_max(j,i) = max(data(td).ORN(whiff_ons(i):whiff_offs(i),j));
		end
		clear i

	end
	clear j

	for i = 1:length(whiff_ons)
		% normalise
		whiff_stim_max(:,i) = whiff_stim_max(:,i)/(whiff_stim_max(1,i));
		whiff_resp_max(:,i) = whiff_resp_max(:,i)/(whiff_resp_max(1,i));
	end

	figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
	set(gcf,'DefaultAxesColorOrder',jet(length(whiff_ons)))
	subplot(1,2,1), hold on
	plot(whiff_stim_max)
	xlabel('Trial Number')
	ylabel('Whiff Peak (norm)')
	set(gca,'YLim',[0 2])

	subplot(1,2,2), hold on
	plot(whiff_resp_max)
	xlabel('Trial Number')
	ylabel('Max Whiff Response (norm)')
	set(gca,'YLim',[0 2])

	PrettyFig;
	snapnow;
	delete(gcf);

end


%%
% This captures what we see before, that even though the stimulus drops precipitously from trial to trial, the neuron response seems relatively unchanged.



% figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% for i = 1:length(whiff_ons)
% 	% calculate stim decay slopes
% 	f1 = fit((1:width(data(td).PID))',whiff_stim_max(:,i),'poly1');
% 	f2 = fit((1:width(data(td).ORN))',whiff_resp_max(:,i),'poly1');
% 	p1 = f1.p1/max(whiff_stim_max(:,i));
% 	p2 = f2.p1/max(whiff_resp_max(:,i));
% 	plot(p1,p2,'.k','MarkerSize',24)

% end


%        ######  ########    ###    ######## ####  ######  ######## ####  ######   ######  
%       ##    ##    ##      ## ##      ##     ##  ##    ##    ##     ##  ##    ## ##    ## 
%       ##          ##     ##   ##     ##     ##  ##          ##     ##  ##       ##       
%        ######     ##    ##     ##    ##     ##   ######     ##     ##  ##        ######  
%             ##    ##    #########    ##     ##        ##    ##     ##  ##             ## 
%       ##    ##    ##    ##     ##    ##     ##  ##    ##    ##     ##  ##    ## ##    ## 
%        ######     ##    ##     ##    ##    ####  ######     ##    ####  ######   ######  

%% Data Statistics and Linear Fit
% The following figure describes the statistics of the stimulus and the response. Left panel: Histograms of stimulus and response. Middle panel: Autocorrelation functions of the stimulus and the response. Right: Linear filter extracted from this dataset. 

% build a simple linear model
if isvector(data(td).PID)
	[K,~,filtertime] = FindBestFilter(data(td).PID,data(td).ORN,[],'filter_length=201;','min_cutoff = 0;');
else
	pid = mean(data(td).PID,2);
	orn = mean(data(td).ORN,2);
	[K,~,filtertime] = FindBestFilter(pid,orn,[],'filter_length=201;','min_cutoff = 0;','offset=40;');
end
data(td).K = K;
data(td).filtertime = filtertime*mean(diff(data(td).time));
if isvector(data(td).PID)
	data(td).LinearFit = mean(data(td).ORN) + convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);
else
	pid = mean(data(td).PID,2);
	orn = mean(data(td).ORN,2);
	data(td).LinearFit = mean(orn) + convolve(data(td).time,pid,data(td).K,data(td).filtertime);
end
data(td).LinearFit(data(td).LinearFit < 0) = 0;

clear ph
figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
ph(1)=subplot(1,3,1); hold on
ph(2)=subplot(1,3,2); hold on
[act] = PlotDataStatistics(data,td,ph);

set(ph(1),'XScale','log','YScale','log')
set(ph(1),'XMinorTick','on','YTick',[0.01 0.1 0.5 1])

subplot(1,3,3), hold on;
plot(data(td).filtertime,data(td).K,'k','LineWidth',2)
set(gca,'XLim',[min(data(td).filtertime) max(data(td).filtertime)],'box','on')
xlabel('Filter Lag (s)')
title('Filter')
PrettyFig;

snapnow;
delete(gcf);


%     ##       #### ##    ## ########    ###    ########         ######## #### ######## 
%     ##        ##  ###   ## ##         ## ##   ##     ##        ##        ##     ##    
%     ##        ##  ####  ## ##        ##   ##  ##     ##        ##        ##     ##    
%     ##        ##  ## ## ## ######   ##     ## ########         ######    ##     ##    
%     ##        ##  ##  #### ##       ######### ##   ##          ##        ##     ##    
%     ##        ##  ##   ### ##       ##     ## ##    ##         ##        ##     ##    
%     ######## #### ##    ## ######## ##     ## ##     ##        ##       ####    ##    


%%
% The weird shape of the linear filter doesn't bode well for the quality of the linear prediction. The following figure shows the linear fit compared to the data. 

figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,2,1), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[25 35])
set(gca,'YLim',[-0.1 1.1*max(max(data(td).PID))])

subplot(2,2,2), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[38 48])
set(gca,'YLim',[-0.1 1.1*max(max(data(td).PID))])

subplot(2,2,3), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,data(td).LinearFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[25 35])

subplot(2,2,4), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,data(td).LinearFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[38 48])
legend Data LinearFit
PrettyFig;


snapnow;
delete(gcf);


%%
% The rsquare of the linear fit is:

if isvector(data(td).ORN)
	disp(rsquare(data(td).LinearFit,data(td).ORN))
else
	disp(rsquare(data(td).LinearFit,mean2(data(td).ORN)))
end

%        ######      ###    #### ##    ##          ###    ##    ##    ###    ##       
%       ##    ##    ## ##    ##  ###   ##         ## ##   ###   ##   ## ##   ##       
%       ##         ##   ##   ##  ####  ##        ##   ##  ####  ##  ##   ##  ##       
%       ##   #### ##     ##  ##  ## ## ##       ##     ## ## ## ## ##     ## ##       
%       ##    ##  #########  ##  ##  ####       ######### ##  #### ######### ##       
%       ##    ##  ##     ##  ##  ##   ###       ##     ## ##   ### ##     ## ##       
%        ######   ##     ## #### ##    ##       ##     ## ##    ## ##     ## ######## 


%% Gain Analysis: Comparison to Linear Model
% Even though the linear model does a pretty bad job estimating the response, can we compare the response to the linear model output to see if there is a systematic variation of gain? 

figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
clear ph
ph(3)=subplot(1,2,1); hold on; 	axis square
ph(4)=subplot(1,2,2); hold on;	axis square
s = 1; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
step_size = 2*act/10;
if step_size < 0.03
	step_size= 0.03;
else
	step_size = 0.03*floor(step_size/0.03);
end
history_lengths = 0:step_size:2*act;
example_history_length = history_lengths(3);

clear x
if isvector(data(td).ORN)
	x.response = data(td).ORN(s:z);
else
	x.response = mean(data(td).ORN(s:z,:),2);
end
x.prediction = data(td).LinearFit(s:z);
if isvector(data(td).PID)
	x.stimulus = data(td).PID(s:z);
else
	x.stimulus = mean(data(td).PID(s:z,:),2);
end
x.time = data(td).time(s:z);
x.filter_length = 201;

redo_bootstrap = 0;
if redo_bootstrap
	data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
end
clear x
set(ph(4),'YScale','log')


snapnow;
delete(gcf);


%%
% Because the fit is so poor, it's hard to fit lines to these clouds of points. We see the general effect where responses to times where the stimulus is low (green) are systematically under-estimated by the linear model, and the responses to times where the stimulus is high (red) are systematically over-estimated by the linear model. However, perhaps this an artefact of the very poor fit?

%%
% In the following analysis we redo the gain analysis, but restrict the analysis only to the whiffs. 

if whiff_anal

	figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
	clear ph
	ph(3)=subplot(1,2,1); hold on; 	axis square
	ph(4)=subplot(1,2,2); hold on;	axis square
	s = 1; % when we start for the gain analysis
	z = length(data(td).ORN); % where we end
	step_size = 2*act/10;
	if step_size < 0.03
		step_size= 0.03;
	else
		step_size = 0.03*floor(step_size/0.03);
	end
	history_lengths = [0:step_size:2*act];
	example_history_length = history_lengths(3);

	clear x
	if isvector(data(td).ORN)
		x.response = data(td).ORN(s:z);
	else
		x.response = mean(data(td).ORN(s:z,:),2);
	end
	x.prediction = data(td).LinearFit(s:z);
	if isvector(data(td).PID)
		x.stimulus = data(td).PID(s:z);
	else
		x.stimulus = mean(data(td).PID(s:z,:),2);
	end

	% restrict to whiffs
	x.stimulus(~hs) = NaN;

	x.time = data(td).time(s:z);
	x.filter_length = 201;

	redo_bootstrap = 0;
	if redo_bootstrap
		data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
	else
		GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
	end
	clear x



	snapnow;
	delete(gcf);

end


%      ##    ## ##       ##    ##       ##     ##  #######  ########  ######## ##       
%      ###   ## ##       ###   ##       ###   ### ##     ## ##     ## ##       ##       
%      ####  ## ##       ####  ##       #### #### ##     ## ##     ## ##       ##       
%      ## ## ## ##       ## ## ##       ## ### ## ##     ## ##     ## ######   ##       
%      ##  #### ##       ##  ####       ##     ## ##     ## ##     ## ##       ##       
%      ##   ### ##       ##   ###       ##     ## ##     ## ##     ## ##       ##       
%      ##    ## ######## ##    ##       ##     ##  #######  ########  ######## ######## 


%% Fitting a NLN model
% In this section we attempt to explain the data by fitting a Nonlinear-Linear-Nonlinear model. This model consists of a static input non-linearity, a parametric bi-lobed filter, and a static output non-linearity. The static non-linearities are Hill functions.  

%%
% The following figure shows the components of the model.

figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
subplot(1,3,1), hold on
x = 0:1e-3:100;
y = hill(data(td).NLNFitParam(1:3),x);
plot(x,y)
set(gca,'XLim',[0 max(max(data(td).PID))])
xlabel('Stimulus (V)')
ylabel('Function Output (a.u.)')
title('Input Nonlinearity')

subplot(1,3,2), hold on
t = 1:50;
x = data(td).NLNFitParam;
y = make_bilobe_filter(x(4),x(5),x(6),x(7),t);
plot(t*3e-3,y)
xlabel('Lag (s)')
ylabel('Filter amplitude (a.u.)')
title('Filter')

subplot(1,3,3), hold on
x = 0:1e-3:100;
y = hill(data(td).NLNFitParam(8:10),x);
plot(x,y)
xlabel('Input (a.u.)')
ylabel('Function Output (Hz)')
title('Output Nonlinearity')

clear x y t

PrettyFig;


snapnow;
delete(gcf);


%%
% The following figure compares the neuron's response (black) to the NLN prediction (red). 
if isvector(data(td).PID)
	NLNFit = SolveNLNModel(data(td).NLNFitParam,data(td).PID);
else
	NLNFit = SolveNLNModel(data(td).NLNFitParam,mean2(data(td).PID));
end


figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,2,1), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[25 35])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

subplot(2,2,2), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[38 48])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

subplot(2,2,3), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,NLNFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[25 35])

subplot(2,2,4), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,NLNFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[38 48])
legend Data NLNFit
PrettyFig;


%%
% The r-square of the NLN model fit is:

if isvector(data(td).ORN)
	disp(rsquare(NLNFit,data(td).ORN))
else
	disp(rsquare(NLNFit,mean2(data(td).ORN)))
end


figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
clear ph
ph(3)=subplot(1,2,1); hold on; 	axis square
ph(4)=subplot(1,2,2); hold on;	axis square
s = 1; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
step_size = 2*act/10;
if step_size < 0.03
	step_size= 0.03;
else
	step_size = 0.03*floor(step_size/0.03);
end
history_lengths = 0:step_size:2*act;
example_history_length = history_lengths(3);

clear x
if isvector(data(td).ORN)
	x.response = data(td).ORN(s:z);
else
	x.response = mean(data(td).ORN(s:z,:),2);
end
x.prediction = NLNFit(s:z);
if isvector(data(td).PID)
	x.stimulus = data(td).PID(s:z);
else
	x.stimulus = mean(data(td).PID(s:z,:),2);
end
x.time = data(td).time(s:z);
x.filter_length = 201;

redo_bootstrap = 0;
if redo_bootstrap
	data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
end
clear x
set(ph(4),'YScale','log')


snapnow;
delete(gcf);

%%
% Again, we repeat the gain analysis but restrict the analysis only to the whiffs, ignoring the blanks inbetween. 

if whiff_anal

	figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
	clear ph
	ph(3)=subplot(1,2,1); hold on; 	axis square
	ph(4)=subplot(1,2,2); hold on;	axis square
	s = 1; % when we start for the gain analysis
	z = length(data(td).ORN); % where we end
	step_size = 2*act/10;
	if step_size < 0.03
		step_size= 0.03;
	else
		step_size = 0.03*floor(step_size/0.03);
	end
	history_lengths = [0:step_size:2*act];
	example_history_length = history_lengths(3);

	clear x
	if isvector(data(td).ORN)
		x.response = data(td).ORN(s:z);
	else
		x.response = mean(data(td).ORN(s:z,:),2);
	end
	x.prediction = NLNFit(s:z);
	if isvector(data(td).PID)
		x.stimulus = data(td).PID(s:z);
	else
		x.stimulus = mean(data(td).PID(s:z,:),2);
	end

	% restrict to whiffs
	x.stimulus(~hs) = NaN;

	x.time = data(td).time(s:z);
	x.filter_length = 201;

	redo_bootstrap = 0;
	if redo_bootstrap
		data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
	else
		GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
	end
	clear x

end


%      ########     ###             ##     ##  #######  ########  ######## ##       
%      ##     ##   ## ##            ###   ### ##     ## ##     ## ##       ##       
%      ##     ##  ##   ##           #### #### ##     ## ##     ## ##       ##       
%      ##     ## ##     ##          ## ### ## ##     ## ##     ## ######   ##       
%      ##     ## #########          ##     ## ##     ## ##     ## ##       ##       
%      ##     ## ##     ##          ##     ## ##     ## ##     ## ##       ##       
%      ########  ##     ##          ##     ##  #######  ########  ######## ######## 



%% Fitting a DA Model
% Can a DA Model explain the responses of this neuron in this dataset? The following figure shows the ORN firing rates and the best-fit DA Model. 

if isvector(data(td).PID)
	data(td).DAFit=DA_integrate2(data(td).PID,data(td).DAFitParam);
else
	data(td).DAFit=DA_integrate2(mean2(data(td).PID),data(td).DAFitParam);
end


figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,2,1), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[25 35])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

subplot(2,2,2), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',[38 48])
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

subplot(2,2,3), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,data(td).DAFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[25 35])

subplot(2,2,4), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,data(td).DAFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',[38 48])
legend Data DAFit
PrettyFig;



snapnow;
delete(gcf);


%%
% The r-square of the DA model fit is:

if isvector(data(td).ORN)
	disp(rsquare(data(td).DAFit,data(td).ORN))
else
	disp(rsquare(data(td).DAFit,mean2(data(td).ORN)))
end

%%
% and the parameters of the model are:

disp(data(td).DAFitParam)


%% Gain Analysis: Comparison to DA Model
% Does the DA model also show systematic variation of gain for high and low stimuli? 

figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
clear ph
ph(3)=subplot(1,2,1); hold on; 	axis square
ph(4)=subplot(1,2,2); hold on;	axis square
s = 1; % when we start for the gain analysis
z = length(data(td).ORN); % where we end
step_size = 2*act/10;
if step_size < 0.03
	step_size= 0.03;
else
	step_size = 0.03*floor(step_size/0.03);
end
history_lengths = 0:step_size:2*act;
example_history_length = history_lengths(3);

clear x
if isvector(data(td).ORN)
	x.response = data(td).ORN(s:z);
else
	x.response = mean(data(td).ORN(s:z,:),2);
end
x.prediction = data(td).DAFit(s:z);
if isvector(data(td).PID)
	x.stimulus = data(td).PID(s:z);
else
	x.stimulus = mean(data(td).PID(s:z,:),2);
end
x.time = data(td).time(s:z);
x.filter_length = 201;

redo_bootstrap = 0;
if redo_bootstrap
	data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
end
clear x


snapnow;
delete(gcf);

%%
% Again, we repeat the gain analysis but restrict the analysis only to the whiffs, ignoring the blanks inbetween. 


if whiff_anal

	figure('outerposition',[0 0 1000 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
	clear ph
	ph(3)=subplot(1,2,1); hold on; 	axis square
	ph(4)=subplot(1,2,2); hold on;	axis square
	s = 1; % when we start for the gain analysis
	z = length(data(td).ORN); % where we end
	step_size = 2*act/10;
	if step_size < 0.03
		step_size= 0.03;
	else
		step_size = 0.03*floor(step_size/0.03);
	end
	history_lengths = [0:step_size:2*act];
	example_history_length = history_lengths(3);

	clear x
	if isvector(data(td).ORN)
		x.response = data(td).ORN(s:z);
	else
		x.response = mean(data(td).ORN(s:z,:),2);
	end
	x.prediction = data(td).DAFit(s:z);
	if isvector(data(td).PID)
		x.stimulus = data(td).PID(s:z);
	else
		x.stimulus = mean(data(td).PID(s:z,:),2);
	end

	% restrict to whiffs
	x.stimulus(~hs) = NaN;

	x.time = data(td).time(s:z);
	x.filter_length = 201;

	redo_bootstrap = 0;
	if redo_bootstrap
		data(td).LinearFit_p = GainAnalysis3(x,history_lengths,example_history_length,ph);
	else
		GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
	end
	clear x



	snapnow;
	delete(gcf);



end

	