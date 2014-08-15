% Mahmut_Data_Analyis.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% specify which data to look at
clearvars -except options
load('/local-data/DA-paper/mahmut_data.mat')
td = 5;
redo_bootstrap = 0;
whiff_anal = 0;

roi1 = [20 30];
roi2 = [55 65];

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

% assuming baseline has already been removed...


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
a = max([a 1e-6]);
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
set(gca,'XLim',roi1)
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

subplot(3,2,2), hold on
plot(data(td).time,data(td).PID);
ylabel('PID (V)')
set(gca,'XLim',roi2)
set(gca,'YLim',[0 1.1*max(max(data(td).PID))])

if isfield(data,'spiketimes') && isfield(data,'spiketimeB')
	subplot(3,2,3), hold on
	raster2(data(td).spiketimes)
	set(gca,'XLim',roi1)


	subplot(3,2,4), hold on
	raster2(data(td).spiketimes)
	set(gca,'XLim',roi2)
end

subplot(3,2,5), hold on
plot(data(td).time,data(td).ORN);
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',roi1)
set(gca,'YLim',[0 1.1*max(max(data(td).ORN))])

subplot(3,2,6), hold on
plot(data(td).time,data(td).ORN);
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',roi2)
set(gca,'YLim',[0 1.1*max(max(data(td).ORN))])
PrettyFig;

snapnow;
delete(gcf);

%%
% It looks like the amplitude of the signal is dropping every trial. It's still not clear why this is the case. 




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
	[K,~,filtertime] = FindBestFilter(pid,orn,[],'filter_length=201;','min_cutoff = 0;');
end
data(td).K = K;
data(td).filtertime = filtertime*mean(diff(data(td).time));
if isvector(data(td).PID)
	data(td).LinearFit = 0*mean(data(td).ORN) + convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);
else
	pid = mean(data(td).PID,2);
	orn = mean(data(td).ORN,2);
	data(td).LinearFit = 0*mean(orn) + convolve(data(td).time,pid,data(td).K,data(td).filtertime);
end


clear ph
figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
ph(1)=subplot(1,3,1); hold on
ph(2)=subplot(1,3,2); hold on
act = PlotDataStatistics(data,td,ph);

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


%% Linear Fit
% The linear filter calculated and shown above is a "best" filter that simultaneously tries to maximise the predictive power of the model while suppressing high frequency components and having a unit gain. 

% Is the linear filter being calculated correctly? To check, we pass the stimulus through a known filter and then back out the filter from the synthetic data and compare the estimated filter to the known filter. 

% make a fake filter
t = 1:201;
fakeK = 200*make_bilobe_filter(10,2,20,2,t);

% generate fake response
if isvector(data(td).PID)
	PID = data(td).PID;
	ORN = data(td).ORN;
else
	PID = mean2(data(td).PID);
	ORN = mean2(data(td).ORN);
end


fakeORN = filter(fakeK,1,PID - mean(PID));
fakeORN(fakeORN<0)=0;

% extract  filter 
[Kprime,~,ft] = FindBestFilter(PID,fakeORN,[],'filter_length=201;');

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(t,fakeK)
plot(ft,Kprime,'r')
legend KnownFilter ReconstructedFilter
set(gca,'XLim',[min(filtertime) max(filtertime)])
xlabel('Lag (a.u.)')
ylabel('Filter height (a.u.)')

PrettyFig;
snapnow;
delete(gcf);

%%
% The reconstructed filter is approximately correct but not exact. Note that both the temporal shape and the raw amplitude of the known filter are backed out correctly by the filter estimation. Why is it not exact? Probably because the stimulus is sparse. Here, we back out the known filter using random noise as an input. 

PID = (randn(length(PID),1));
fakeORN = filter(fakeK,1,PID - mean(PID));


% extract  filter 
[Kprime,~,ft] = FindBestFilter(PID,fakeORN,[],'filter_length=201;');

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(t,fakeK)
plot(ft,Kprime,'r')
legend KnownFilter ReconstructedFilter
set(gca,'XLim',[min(filtertime) max(filtertime)])
xlabel('Lag (a.u.)')
ylabel('Filter height (a.u.)')


PrettyFig;
snapnow;
delete(gcf);

%%
% So now we can exactly match both the amplitude and the temporal structure of the known filter. So our filter estimation techniques work well. The reason the filter is so bad in the real data is probably because of very sparse data, and the fact that the ORN response can't fully be explained by a linear filter. 


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
set(gca,'XLim',roi1)
set(gca,'YLim',[-0.1 1.1*max(max(data(td).PID))])

subplot(2,2,2), hold on
if isvector(data(td).PID)
	plot(data(td).time,data(td).PID,'k');
else
	plot(data(td).time,mean2(data(td).PID),'k');
end
ylabel('PID (V)')
set(gca,'XLim',roi2)
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
set(gca,'XLim',roi1)

subplot(2,2,4), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,data(td).LinearFit,'r')
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
set(gca,'XLim',roi2)
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


% % Perhaps the reason the filter looks weird is because the stimulus actually is dropping from trial to trial, and that averaging the data over trials is not reasonable. In the following section, we combine all the data in one long vector and then calculate the filter from this using a fixed regularisation. 

% s = data(td).PID(:);
% f = data(td).ORN(:);
% [K,d] = FindBestFilter(s,f,[],'filter_length=201;');

% figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% plot(data(td).filtertime,K,'r')
% set(gca,'XLim',[min(data(td).filtertime) max(data(td).filtertime)])
% xlabel('Lag (a.u.)')
% ylabel('Filter height (Hz)')

% fp = convolve(data(td).time,s,K,data(td).filtertime);

% % The r-square of this prediction is:

% disp(rsquare(f,fp))

%              ##    ## ##                  ######## ##     ## ##    ##  ######  
%              ###   ## ##                  ##       ##     ## ###   ## ##    ## 
%              ####  ## ##                  ##       ##     ## ####  ## ##       
%              ## ## ## ##                  ######   ##     ## ## ## ## ##       
%              ##  #### ##                  ##       ##     ## ##  #### ##       
%              ##   ### ##                  ##       ##     ## ##   ### ##    ## 
%              ##    ## ########            ##        #######  ##    ##  ######  



%% Fitting the Data with an Input Nonlinearity
% The filter doesn't seem to have any obvious shape. Can we fit the data with just a nonlinear transformation of the input? 

if isvector(data(td).ORN)
	d.response = data(td).ORN;
else
	d.response = mean2(data(td).ORN);
end

if isvector(data(td).PID)
	d.stimulus = data(td).PID;

else
	d.stimulus = mean2(data(td).PID);
end

d.stimulus = d.stimulus - min(d.stimulus);
NFit = hill(data(td).NFit,d.stimulus);

figure('outerposition',[0 0 1200 500],'PaperUnits','points','PaperSize',[1200 500]); hold on
subplot(1,3,1), hold on
plot(d.stimulus,hill(data(td).NFit,d.stimulus),'k+')
xlabel('Stimulus (V)')
ylabel('Firing rate (Hz)')
title('Nonlinear Function')

subplot(1,3,2:3), hold on
plot(data(td).time,d.response,'k')
hold on
plot(data(td).time,NFit,'r')
legend({'Data','Prediction'})
set(gca,'XLim',[min(data(td).time) max(data(td).time)])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
title('Nonlinear Function Output')

PrettyFig;
snapnow;
delete(gcf);


%%
% The r-square of the prediction is:
disp(rsquare(NFit,d.response))



%           ##    ## ##          ##     ##  #######  ########  ######## ##       
%           ###   ## ##          ###   ### ##     ## ##     ## ##       ##       
%           ####  ## ##          #### #### ##     ## ##     ## ##       ##       
%           ## ## ## ##          ## ### ## ##     ## ##     ## ######   ##       
%           ##  #### ##          ##     ## ##     ## ##     ## ##       ##       
%           ##   ### ##          ##     ## ##     ## ##     ## ##       ##       
%           ##    ## ########    ##     ##  #######  ########  ######## ######## 



%% Fitting Data with a NL model
% Clearly, the linear prediction isn't doing a very good job. Why? Can we improve the prediction by inserting an input non-linearity into the model? The following figure shows the result of fitting an input non-linearity (top left) and a linear filter (top right) to the data. The lower panel shows the response of the neuron compared to this prediction. 


if isvector(data(td).ORN)
	d.response = data(td).ORN;
else
	d.response = mean2(data(td).ORN);
end

if isvector(data(td).PID)
	d.stimulus = data(td).PID;
else
	d.stimulus = mean2(data(td).PID);
end

d.stimulus = abs(d.stimulus);

[NLFit, K] = SolveNLModel(data(td).NLFitParam,d.stimulus,d.response);
figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,2,1), hold on
scatter(d.stimulus,hill2(data(td).NLFitParam,d.stimulus),'k')
xlabel('Stimulus (V)')
ylabel('Function Output (V)')

set(gca,'XLim',[0 max(d.stimulus)],'YLim',[0 max(d.stimulus)])


subplot(2,2,2), hold on
ft = mean(diff(data(td).time))*(1:length(K));
plot(ft,K)
set(gca,'XLim',[0 max(ft)])
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (Hz)')


subplot(2,2,3:4), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,NFit,'g')
plot(data(td).time,NLFit,'b')
set(gca,'XLim',[20 40])
legend({'Data','Nonlinear Fit','NL Fit'},'Location','EastOutside')
xlabel('Time (s)')
ylabel('Firing rate (Hz)')


PrettyFig
snapnow;
delete(gcf);

%%
% The rsquare of this NL fit is:

if isvector(data(td).ORN)
	disp(rsquare(NLFit,data(td).ORN))
else
	disp(rsquare(NLFit,mean2(data(td).ORN)))
end


%               ##       ##    ##        ##     ##  #######  ########  ######## ##       
%               ##       ###   ##        ###   ### ##     ## ##     ## ##       ##       
%               ##       ####  ##        #### #### ##     ## ##     ## ##       ##       
%               ##       ## ## ##        ## ### ## ##     ## ##     ## ######   ##       
%               ##       ##  ####        ##     ## ##     ## ##     ## ##       ##       
%               ##       ##   ###        ##     ## ##     ## ##     ## ##       ##       
%               ######## ##    ##        ##     ##  #######  ########  ######## ######## 




%% LN Model
% Does putting the static non-linearity at the output of a linear filter improve the prediction? In the following section we find a best-fit NL model, where the linear filter and the static non-linearity are fit simultaneously by a iterative fitting procedure. 

%%
% We use a Hill function for the static non-linearity here (as we do so in the previous cases). Because the Hill function is not defined for negative values, this model actually has a rectifier between the linear output and the output non-linearity: only positive values from the linear block are considered as valid inputs to the non-linearity. 

if isvector(data(td).ORN)
	d.response = data(td).ORN;
else
	d.response = mean2(data(td).ORN);
end

if isvector(data(td).PID)
	d.stimulus = data(td).PID;
else
	d.stimulus = mean2(data(td).PID);
end


[LNFit, K] = SolveLNModel(data(td).LNFitParam,d.stimulus,d.response);
figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
subplot(2,2,2), hold on
scatter(ihill2(data(td).LNFitParam,d.response),d.response,'k')
xlabel('Filter Output (a.u.)')
ylabel('Neuron Response (Hz)')

set(gca,'XLim',[0 max(d.response)],'YLim',[0 max(d.response)])


subplot(2,2,1), hold on
ft = mean(diff(data(td).time))*(1:length(K));
plot(ft,K)
set(gca,'XLim',[0 max(ft)])
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (Hz)')


subplot(2,2,3:4), hold on
if isvector(data(td).ORN)
	plot(data(td).time,data(td).ORN,'k');
else
	plot(data(td).time,mean2(data(td).ORN),'k');
end
plot(data(td).time,NFit,'g')
plot(data(td).time,NLFit,'b')
plot(data(td).time,LNFit,'r')
set(gca,'XLim',[20 40])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Data','Nonlinear Fit','NL Fit','LN Fit'},'Location','EastOutside')

PrettyFig
snapnow;
delete(gcf);

%%
% The rsquare of this LN fit is:

if isvector(data(td).ORN)
	disp(rsquare(LNFit,data(td).ORN))
else
	disp(rsquare(LNFit,mean2(data(td).ORN)))
end


return




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



%    ##      ## ##     ## #### ######## ########          ###    ##    ##    ###    ##       
%    ##  ##  ## ##     ##  ##  ##       ##               ## ##   ###   ##   ## ##   ##       
%    ##  ##  ## ##     ##  ##  ##       ##              ##   ##  ####  ##  ##   ##  ##       
%    ##  ##  ## #########  ##  ######   ######         ##     ## ## ## ## ##     ## ##       
%    ##  ##  ## ##     ##  ##  ##       ##             ######### ##  #### ######### ##       
%    ##  ##  ## ##     ##  ##  ##       ##             ##     ## ##   ### ##     ## ##       
%     ###  ###  ##     ## #### ##       ##             ##     ## ##    ## ##     ## ######## 



%% Whiff-based Analysis
% To quantify this, let's attempt to find all the times when the stimulus is high, i.e., ten standard deviations above the baseline. This threshold ensures that we pick up all the whiffs, but ignore the blanks in between, as shown in the figure below. 

if whiff_anal

	figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
	plot(data(td).time, data(td).PID(:,1),'k'), hold on

	hs = data(td).PID(:,1) > 10*std(data(td).PID(1:3000,1));
	plot(data(td).time(hs), data(td).PID(hs,1),'.r','MarkerSize',24), hold on
	set(gca,'XLim',roi1)
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




	