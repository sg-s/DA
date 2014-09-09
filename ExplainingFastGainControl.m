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

% load data
load('/local-data/DA-paper/data.mat')
td = 7;

%% Fitting a DA Model to 1-octen-3-ol data
% In this section, we fit a DA model to the example 1-octen-3-ol data: 

disp(data(td).original_name)

% fit model to data
if ~exist('p')

	% choose model
	global DA_Model_Func
	global DA_Model_Validate_Param_Func
	global nsteps
	nsteps = 80;
	DA_Model_Func = @DA_integrate2;
	DA_Model_Validate_Param_Func = @ValidateDAParameters2;

	clear d
	d.stimulus = data(td).PID/max(data(td).PID);
	d.response = data(td).ORN/max(data(td).ORN);
	[p,~,x] = FitDAModelToData(d,[5.2 6.6 0 .78 11 6.4 6 .557],[1 1 0 1e-1 1 0 1 -1],[10 20 1 20 20 24 10 1]);

	DApred=DA_Model_Func(d.stimulus,p)*max(data(td).ORN);
	multiplot(data(td).time,data(td).ORN,DApred)
end

DApred=DA_integrate2(data(td).PID,p);

% fit a LN model (re-use old fits)
load('HowGeneralIsGainAdaptation.mat','Filters','HillFit','filtertime','history_lengths','p_values_low')
LinearFit = mean(data(td).ORN)+convolve(data(td).time,data(td).PID,Filters(:,td),filtertime);
LinearFit(LinearFit<0)=0;
LNpred = hill(HillFit(:,td),LinearFit);

% make a plot showing the data, the LN fit, and the DA fit
fh=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID,'k');
set(gca,'XLim',[mean(data(td).time)-2 mean(data(td).time)+2])
ylabel('PID Voltage (V)')
hold off


subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(td).time,LNpred,'g');
plot(data(td).time,DApred,'r')
set(gca,'XLim',[mean(data(td).time)-2 mean(data(td).time)+2])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
legend({'Data','LN Fit','DA Fit'})
PrettyFig;
hold off


if being_published
	snapnow;
	delete(fh);
end

%%
% The r-square of the LN model prediction is: 

s = 300;
z = length(data(td).ORN)-33;
disp(rsquare(data(td).ORN(s:z),LNpred(s:z)))

%%
% _cf._ r-square of the DA model prediction: 

s = 300;
z = length(data(td).ORN)-33;
disp(rsquare(data(td).ORN(s:z),DApred(s:z)))

%%
% The following figure shows the best-fit LN model. The panel on the left shows the linear filter, and the panel on the right shows the output non-linearity. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(filtertime,Filters(:,td),'k')
xlabel('Lag (s)')
ylabel('Filter Amplitude (Hz)')
set(gca,'XLim',[min(filtertime) max(filtertime)])
subplot(1,2,2), hold on
x = 0:1:max(data(td).ORN);
y = hill(HillFit(:,td),x);
plot(x,y,'k')
xlabel('Filter Output (Hz)')
ylabel('Nonlinearity Output (Hz)')
PrettyFig;

if being_published
	snapnow;
	delete(fh);
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
	p_DA=GainAnalysis3(x,history_lengths,example_history_length,ph);
else
	p_DA=GainAnalysis3(x,history_lengths,example_history_length,ph,p_DA);
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