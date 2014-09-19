% ExplainingFastGainControl.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Explaining Fast Gain Control
% We have shown that there is a lot of evidence that ORNs modulate gain as a function of recnetly experienced stimuli. We show that the LN model, despite explaining >95% of the data (typically), fails to account for this fast gain control. Can we develop a model to explain this fast gain control?

%%
% In the following sections, we try fitting a Dynamical Adaptation (DA) model to the data. 

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end


% some parameters
font_size = 20;
marker_size = 10;
marker_size2 = 20;
history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;
example_history_length = 0.135;


% load data
load('/local-data/DA-paper/data.mat')
td = 4;


%         #######  ##     ## ######## ########  ##     ## ########    ##    ## ##       
%        ##     ## ##     ##    ##    ##     ## ##     ##    ##       ###   ## ##       
%        ##     ## ##     ##    ##    ##     ## ##     ##    ##       ####  ## ##       
%        ##     ## ##     ##    ##    ########  ##     ##    ##       ## ## ## ##       
%        ##     ## ##     ##    ##    ##        ##     ##    ##       ##  #### ##       
%        ##     ## ##     ##    ##    ##        ##     ##    ##       ##   ### ##       
%         #######   #######     ##    ##         #######     ##       ##    ## ######## 

%% Output Non-linearity
% Does adding an output non-linearity post-hoc improve the linear fit? In the following section, we fit a two-parameter Hill function to the linear prediction and the data, after we have generated the linear prediction. (In other words, the Hill function is fit _after_ the best linear prediction is calculated.) The following figure shows the filter and the shape of the best-fit non-linearity:


% build a simple linear model
[K,~,filtertime] = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;');
data(td).K = K;
data(td).filtertime = filtertime*mean(diff(data(td).time));
data(td).LinearFit = convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);
data(td).LinearFit = data(td).LinearFit + mean(data(td).ORN);

xdata = data(td).LinearFit;
ydata = data(td).ORN;

% crop it to lose NaNs
ydata(isnan(xdata)) = [];
xdata(isnan(xdata)) = [];

xdata = xdata(:);
ydata = ydata(:);

fo=optimset('MaxFunEvals',1000,'Display','none');
x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);
% save this for later
LNpred = hill(x,data(td).LinearFit);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(data(td).filtertime,K,'k')
set(gca,'XLim',[min(data(td).filtertime)-.1 max(data(td).filtertime)+.1])
xlabel('Filter Lag (s)')
ylabel('Filter amplitude (Hz/stim)')
subplot(1,2,2), hold on
plot(sort(xdata),hill(x,sort(xdata)),'k')
xlabel('Linear Prediction (Hz)')
ylabel('Nonlinearity Output (Hz)')
set(gca,'XLim',[0 max(xdata)])

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end


%%
% How does the output nonlinearity change the prediciton? In the following figure, the data is shown in black,and the LN prediction is shown in red. 


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,LNpred,'r')
set(gca,'XLim',[18 22])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end




%%
% They look almost identical. The r-square of the LN prediction is: 

disp(rsquare(LNpred,data(td).ORN))

%%
% And the raw Euclidean Distance between the prediction and the data is: 

disp(Cost2(LNpred(205:end-33),data(td).ORN(205:end-33)))

%      ########     ###       ##     ##  #######  ########  ######## ##       
%      ##     ##   ## ##      ###   ### ##     ## ##     ## ##       ##       
%      ##     ##  ##   ##     #### #### ##     ## ##     ## ##       ##       
%      ##     ## ##     ##    ## ### ## ##     ## ##     ## ######   ##       
%      ##     ## #########    ##     ## ##     ## ##     ## ##       ##       
%      ##     ## ##     ##    ##     ## ##     ## ##     ## ##       ##       
%      ########  ##     ##    ##     ##  #######  ########  ######## ######## 

%% Fitting a DA Model to ORN response data
% In this section, we fit a DA model to the example data: 
disp(data(td).original_name)

% % fit model to data
% if ~exist('p')
% 	clear d
% 	d.stimulus = data(td).PID - mean(data(td).PID);
% 	d.stimulus = d.stimulus*100;
% 	d.response = data(td).ORN;
% 	x0 = [549 4.8 2e-4 .4 14 8 3 -.1127];
% 	[p,DApred,DAParam]=FitDAModelToData(d,x0,[],[],.95);
% 	close all
% 	multiplot(data(td).time,data(td).ORN,DApred)
% end



% fit model to data
if ~exist('p')
	clear d
	d.stimulus = data(td).PID - mean(data(td).PID);
	d.stimulus = d.stimulus*100;
	d.stimulus = d.stimulus(200:end-200);
	d.LNpred = LNpred(200:end-200);
	d.response = data(td).ORN(200:end-200);
	x0 = [1 .32 .2 .4 14 8 3 -.1127];
	[p,DApred,DAParam]=FitDAModelToData(d,x0,[],[],.95);


	% fix DApred
	temp = NaN*LNpred;
	temp(200:end-200) = DApred;
	DApred = temp; clear temp


end

%%
% The DA model primarily consists of two filters, $K_{z}$ and $K_{y}$, which are shown in the figure below: 

t = 0:0.1:100;
[Ky,Kz] = DA_Filters(p,t);
t = 3e-3*t;
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(t,Ky,'k')
plot(t,Kz,'r')
legend({'K_y','K_z'})
xlabel('Lag (s)')
ylabel('Filter Height')

PrettyFig;

if being_published
	snapnow;
	delete(fh);
end

%%
% All the parameters of the best-fit DA model are:
disp(p)


%%
% How does the prediction of the DA model compare with the LN model? The following figure shows the data, with the two predictions:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data(td).time,data(td).ORN,'k')
plot(data(td).time,LNpred,'r')
plot(data(td).time,DApred,'g')
set(gca,'XLim',[18 22])
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
legend({'Data','LN Model','DA Model'})
PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%%
% They look almost identical. The r-square of the DA prediction is: 

disp(rsquare(DApred,data(td).ORN))

%%
% And the raw Euclidean Distance between the prediction and the data is: 

disp(Cost2(DApred(205:end-33),data(td).ORN(205:end-33)))

%  ######      ###    #### ##    ##       ###    ##    ##    ###    ##          ##       ##    ## 
% ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##          ##       ###   ## 
% ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##          ##       ####  ## 
% ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##          ##       ## ## ## 
% ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##          ##       ##  #### 
% ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##          ##       ##   ### 
%  ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ########    ######## ##    ## 


%% LN Model Gain Analysis
% In the previous section, we saw that the simple linear model fails to account for this fast adaptation of the ORNs. Specifically, the gain of the neuron _w.r.t_ to the model is significantly different for times when the stimulus is high and the when the stimulus is low.

%%
% In this section, we want to know if the LN model (adding a non-linear function post-hoc) corrects for this mis-prediction of gain. Here, we repeat the gain analysis as in the previous section, but this time, using the LN prediction instead of the linear prediction. 

f1=figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
ph(1) = subplot(2,1,1); hold on 
ph(2) = subplot(2,1,2); hold on

title(ph(1),strrep(data(td).original_name,'_','-'),'FontSize',20);

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on



s = 300; % when we start for the gain analysis
z = length(data(td).ORN) - 33; % where we end

clear x
x.response = data(td).ORN(s:z);
x.prediction = LNpred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;

if redo_bootstrap
	[p_LN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	s=abs(l-h);
	s(p_LN(1,:)>0.05)=NaN;
	[~,loc]=max(s);
	example_history_length_LN = history_lengths(loc);

else
	GainAnalysis4(x,history_lengths,example_history_length_LN,ph,p_LN);
end

xlabel(ph(3),'LN Prediction (Hz)')
set(ph(4),'XScale','log')

if being_published
	snapnow;
	delete(f1);

	snapnow;
	delete(f2);
end


 % ######      ###    #### ##    ##       ###    ##    ##    ###    ##          ########     ###    
% ##    ##    ## ##    ##  ###   ##      ## ##   ###   ##   ## ##   ##          ##     ##   ## ##   
% ##         ##   ##   ##  ####  ##     ##   ##  ####  ##  ##   ##  ##          ##     ##  ##   ##  
% ##   #### ##     ##  ##  ## ## ##    ##     ## ## ## ## ##     ## ##          ##     ## ##     ## 
% ##    ##  #########  ##  ##  ####    ######### ##  #### ######### ##          ##     ## ######### 
% ##    ##  ##     ##  ##  ##   ###    ##     ## ##   ### ##     ## ##          ##     ## ##     ## 
 % ######   ##     ## #### ##    ##    ##     ## ##    ## ##     ## ########    ########  ##     ## 


%% DA Model Gain Analysis
%  Here, we repeat the gain analysis as in the previous section, but this time, using the DA prediction instead of the LN prediction. 

f1=figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
ph(1) = subplot(2,1,1); hold on 
ph(2) = subplot(2,1,2); hold on

title(ph(1),strrep(data(td).original_name,'_','-'),'FontSize',20);

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on



s = 300; % when we start for the gain analysis
z = length(data(td).ORN) - 200; % where we end

clear x
x.response = data(td).ORN(s:z);
x.prediction = DApred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;

if redo_bootstrap
	[p_DA,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);

else
	GainAnalysis4(x,history_lengths,example_history_length_LN,ph,p_DA);
end

xlabel(ph(3),'LN Prediction (Hz)')
set(ph(4),'XScale','log')

if being_published
	snapnow;
	delete(f1);

	snapnow;
	delete(f2);
end




return

%%
% What do to the two models predict of the stimulus-dependent instantaneous gain of the neuron? On the left is the prediction of the LN model, on the right is the prediction of the DA model. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
example_history_length = 0.3;
s = 300; % when we start for the gain analysis
z = length(data(td).ORN) - 33; % where we end

ph = []; ph(3) = subplot(1,2,1); hold on
clear x
x.response = data(td).ORN(s:z);
x.prediction = LNpred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;
GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
xlabel('LN Prediction (Hz)')
ylabel('ORN Response (Hz)')

ph = []; ph(3) = subplot(1,2,2); hold on
x.prediction = DApred(s:z);
GainAnalysis3(x,history_lengths,example_history_length,ph,NaN*history_lengths);
xlabel('DA Prediction (Hz)')
ylabel('')



%%
% Are these slopes significant? Over what history lengths are they significantly different? 
figure('outerposition',[0 0 1200 500],'PaperUnits','points','PaperSize',[1200 500]); hold on
example_history_length = 0.12;
s = 300; % when we start for the gain analysis
z = length(data(td).ORN) - 33; % where we end

ph = []; ph(4) = subplot(1,2,1); hold on
clear x
x.response = data(td).ORN(s:z);
x.prediction = LNpred(s:z);
x.stimulus = data(td).PID(s:z);
x.time = data(td).time(s:z);
x.filter_length = 201;
GainAnalysis3(x,history_lengths,example_history_length,ph,p_values(:,td));
set(gca,'XScale','log')
xlabel('History Length (s)')
title('LN Model')


ph = []; ph(4) = subplot(1,2,2); hold on
x.prediction = DApred(s:z);
if ~exist('p_DA')
	p_DA=GainAnalysis4(x,history_lengths,example_history_length,ph);
else
	p_DA=GainAnalysis4(x,history_lengths,example_history_length,ph,p_DA);
end
xlabel('History Length (s)')
ylabel('')
set(gca,'XScale','log')
title('DA Model')

PrettyFig;

if being_published
	snapnow;
	delete(fh);
end


%% The DA Model cannot fit the data if $r_{0}$ >0
% The $r_{0}$ parameter is a fudge factor added to fit some supposed offset which corresponds to the basal firing of the ORN. However, in this case, a negative value is chosen, which makes no sense. If we constrain this to positive values, the best fit is horrible: 

% fit model to data
if ~exist('p2')

	% choose model
	global DA_Model_Func
	global DA_Model_Validate_Param_Func
	global nsteps
	nsteps = 200;
	DA_Model_Func = @DA_integrate2;
	DA_Model_Validate_Param_Func = @ValidateDAParameters2;

	clear d
	d.stimulus = data(td).PID;
	d.response = data(td).ORN;
	[p2,~,x] = FitDAModelToData(d,[7631 105 .003 0.4 7 4.5 10 1],[1 1 0 1e-1 2 0 2 0],[1e5 1e3 1 20 10 24 15 100]);
	clear d
	DApred=DA_Model_Func(data(td).PID,p2);
	multiplot(data(td).time,data(td).ORN,DApred)
end

DApred=DA_integrate2(data(td).PID,p2);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID,'k');
set(gca,'XLim',[min(data(td).time)-2 max(data(td).time)+2])
ylabel('PID Voltage (V)')
hold off


subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(td).time,DApred,'r')
set(gca,'XLim',[min(data(td).time)-2 max(data(td).time)+2])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
legend({'Data','DA Fit'})
PrettyFig;
hold off


if being_published
	snapnow;
	delete(fh);
end


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))