%% Dynamical Adaptation in ORNs 
% Do ORNs exhibit fast adaptation to a flickering stimulus? Can a simple dynamical adaptation model predict ORN responses to flickering inputs? Here, I take data from Carlotta's random stimulus experiments and first check how well a simple linear prediction from the stimulus compares to the data. Then, I study how the instantaneous gain of the actual ORN response w.r.t to the predicted ORN response varies with time, and try to find a objective filter from the stimulus to this instantaneous gain to correct for systematic errors in the linear prediction.  

%%

% parameters to tune for figure display
font_size = 20;
marker_size = 10;
marker_size2 = 20;

%% Data Visualisation 
% Data from this file will be used for this analysis.

filename = '/data/random-stim/final_2011_06_14_ab3A_1o3ol3X-3_20ml_30sec_30ms_rand.mat';
disp(filename)

%%
% Here, data from simultaneous measurements of the stimulus and from the ORN is shown. The top row shows the valve command signal (a binary signal, high means valve is open). The middle row shows the simultaneous measurement of the odour stimulus with a PID, and the bottom row shows the instantaneous firing rate of the ORN recorded from. 

crop_traces = 0;
% grab data
if ~(exist('PID') == 1)
	crop_traces = 1;
	[PID time f fs Valve ntrials] = PrepData2(filename);

	% detrend PID
	ptrend = fit(time',PID','Poly1'); 
	PID = PID - (ptrend(time) - mean(ptrend(time)))';

end

% plot the data
figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(3,1,1), hold on
plot(time,Valve,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[18 22],'YLim',[-0.1 1.1])
ylabel('stimulus')

subplot(3,1,2), hold on
plot(time,PID,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[18 22])
ylabel('PID (a.u.)')

subplot(3,1,3), hold on
plot(time,f,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[18 22])
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

%% 
% The odour used, the neuron recorded from and the correlation time of the flickering stimulus are in the file name displayed above the plot.

%% Filter Extraction 
% Now, we extract a simple linear filter from the stimulus such that a convolution with the input (stimulus) gives the output (ORN firing rate). 

%%
% The filter is calculated using Chichilnisky's method. The only modification I have used is to regularise the matrix before inversion. The code is based on what Carlotta uses, and I have re-written it for clarity and modified the regularisation step.

%% 
% The filter _K_ is given by 
%
% $$ C*K=s'*f $$
%
% where C is 
%
% $$ C=Cov(s)+\frac{rI}{1+r} $$
% 
% where _I_ is the identity matrix and _r_ is a free parameter called the regularisation factor that suppresses the high-frequency components of _K_. _s_ is the stimulus vector (e.g. the PID) and _f_ is the response vector (here, the firing rate of the ORN). In practice, _K_ is estimated by a left matrix division:
%
% $$ K=C\setminus(s'*f) $$

%%
% For comparison, I have also calculated the filter using Damon's function (Back_out_1dfilter_new)

%%
% The problem of choosing _r_ appropriately remains. If _r_ is too low, then the effective number of parameters of the filter goes up, permitting over-fitting of the data, and you end up with a noisy squiggle of a filter that fits the data well. If _r_ is too high, the filter has too few modes, and doesn't fit the data well, but looks nice. The ideal _r_ is chosen so that it is equal to 
%
% $$ \arg\min_{r}\left\Vert \frac{e(r)-\min(e)}{\min(e)},\frac{S(r)-\min(S)}{\min(S)},\frac{\left|m(r)-1\right|}{1}\right\Vert _{\infty} $$
% 
% where _e_ is the error of the prediction w.r.t to the actual data, _S_ is the absolute sum of the filter _K_, and _m_ is the slope of the best fit of the prediction to the data. Minimising _e_ minimises the error of the fit. Minimising _S_ minimises the high-frequency components in the filter, and prevents over-fitting. We also want the slope _m_ to be as close to 1 as possible. 
%%
% For a detailed look at how the slope of best fit and the filter shape vary with _r_, look at the next section.


% shift input, black magic. why are we shifting input? no clue. but this is what Carlotta does. is there a lag between the to traces? why is it this? 
shift_input = 33;
if crop_traces
	time = time(1:end-shift_input+1);
	Valve = Valve(1:end-shift_input+1);
	PID = PID(shift_input:end);
	f = f(1:end-shift_input+1);
	fs = fs(1:end-shift_input+1);
end

% compute filter
filter_length = 333;
filtertime = 0:mean(diff(time)):filter_length*mean(diff(time)); % this is the x axis for the filter K
% both damon's and my filter functions have a free parameter. so we calculate the filter, and find the best value of the free parameter that minimises $ \Sigma |K| $, which will yield the cleanest filters with the least amount of high frequency noise. 
if ~(exist('K') == 1)
	[K Kdamon diagnostics] = FindBestFilter(PID,f,filter_length);
end


% now plot the filter
%% 
% In the figure below, the filter computed the PID and the firing rate is shown. In blue is the filter calculated using Chichilnisky's method, and in red is the filter calculated using Damon's function. The x axis is the lag of the filter, and the Y axis is the amplitude of the filter (in Hz)
figure('outerposition',[0 0 350 350],'PaperUnits','points','PaperSize',[800 350]); hold on
plot(filtertime,K,'LineWidth',2)
plot(filtertime(2:end),Kdamon,'r','LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude (Hz)','FontSize',font_size)
legend Chichilnisky Damon

%%
% The convolution of this filter with the input (PID) gives a prediction of the output. In the figure below, the data is shown in black, and the prediction from this linear filter is shown in blue (Chichilnisky method) and red (Damon's function). Note that sometimes the prediction is below zero, which doesn't mean anything physical (firing rate cannot be < 0). 

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on
fp = filter(K,1,PID);
% censor initial prediction
fp(1:filter_length+1)=NaN;
fpd = filter(Kdamon,1,PID-mean(PID))+mean(f);
% censor initial prediction
fpd(1:filter_length+1)=NaN;
plot(time,f,'k','LineWidth',2)
plot(time,fp,'b','LineWidth',2)
plot(time,fpd,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-2 mean(time)+2],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
legend Data Chichilnisky Damon


%% Analysis of Linear Prediction - Linearity of Prediction 
% For a linear filter calculated from the data, a plot of the actual response to the predicted response can be fit with a line of slope 1. Here the actual firing rate is plotted on the X axis and the firing rate predicted from the linear filter is plotted on the Y axis. 



figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(fp,f,'.','MarkerSize',marker_size)
set(gca,'LineWidth',2,'FontSize',20,'box','on')
set(gca,'XLim',[min([fp f]) max([fp f])],'YLim',[min([fp f]) max([fp f])])
axis square
ylabel('Actual Firing rate (Hz)','FontSize',20)
xlabel('Linear Predictionf Firing rate (Hz)','FontSize',20)
[fall gof] = fit(fp(filter_length+2:end)',f(filter_length+2:end)','Poly1');
plot(fall(min(f):1:max(f)),min(f):1:max(f),'k','LineWidth',2)
title('Chichilnisky','FontSize',font_size)
subplot(1,2,2), hold on
plot(fpd,f,'.','MarkerSize',marker_size)
set(gca,'LineWidth',2,'FontSize',20,'box','on')
set(gca,'XLim',[min([fpd f]) max([fpd f])],'YLim',[min([fpd f]) max([fpd f])])
axis square
ylabel('Actual Firing rate (Hz)','FontSize',20)
xlabel('Linear Prediction Firing rate (Hz)','FontSize',20)
[fall2 gof2] = fit(fpd(filter_length+2:end)',f(filter_length+2:end)','Poly1');
plot(fall2(min(f):1:max(f)),min(f):1:max(f),'k','LineWidth',2)
title('Damon','FontSize',font_size)

%% 
% The line is a best fit to all the data. The slope of the best fit line to Chichilnisky's fit is: 
disp(fall.p1)

%% 
% The line is a best fit to all the data. The slope of the best fit line to Damon's fit is: 
disp(fall2.p1)


%%
% Note that for the Chichilnisky filter, the slope is very close to 1. For very little regularisation, the slope is exactly 1.

%%
% Both Damon's and my filter estimation functions have a free parameter, _r_, the regularisation factor. The following plots show how the error, the filter sum, and the filter height vary with varying the free parameter _r_. The point in red indicates the value of _r_ finally chosen for the best filter. 
PlotFilterDiagnostics(diagnostics);


%% Analysis of Linear Prediction - Which input predicts output better?
% The measurement of the stimulus should help us predict the response of the neuron better than just the valve control signal, as we now capture the fast fluctuations of the odour stimulus as it reaches the ORN. Thus, a prediction of ORN firing response from the PID signal should be better than a prediction of the ORN firing response from the valve signal. 

%% 
% First, we recalculate a new filter from the valve control signal and the ORN firing rate using the same methods as before. The plot below shows this new filter:

figure('outerposition',[0 0 350 350],'PaperUnits','points','PaperSize',[800 350]); hold on
if ~(exist('K2') == 1)
	[K2 K2damon diagnostics_valve] = FindBestFilter(Valve',f,filter_length);
end
plot(filtertime,K2,'LineWidth',2), hold on
plot(filtertime(2:end),K2damon,'r','LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude (Hz)','FontSize',font_size)
title('Filter: Valve->ORN','FontSize',20)
legend Chichilnisky Damon

%%
% Using the filter calculated using Damon's code, we make a new prediction of the firing rate of the ORN. In the figure below, the data is shown in blue, the linear prediction from the PID is shown again in red, and the linear prediction from the Valve is shown in black. 

figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
fp2 = filter(K2damon,1,Valve-mean(Valve)) + mean(f);
fp2(1:filter_length+1) = NaN;
plot(time,f,'LineWidth',2)
plot(time,fp,'r','LineWidth',2)
plot(time,fp2,'k','LineWidth',2)
set(gca,'XLim',[mean(time)-2 mean(time)+3],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
legend Data PIDPrediction ValvePrediction 

%% 
% The quality of the prediction is measured using the Root-mean-square error, calculated as:
%
% $$ \triangle t*\sqrt[]{\Sigma \left(f-f_{p}\right)^{2}} $$
% 
% The error in the prediction from the PID, per unit time is: 

disp(sqrt(sum( (fp(filter_length+2:end) - f(filter_length+2:end)).^2 ))*mean(diff(time)));

%%
% while the error for the prediction from the Valve is: 
fp2 = fp2';
disp(sqrt(sum( (fp2(filter_length+2:end) - f(filter_length+2:end)).^2 ))*mean(diff(time)));


%% Analysis of Linear Prediction - Response to High and Low Stimuli
% How does the response of the ORN differ for high and low stimuli? Specifically, does the neuron display the characteristics of fast adaptation to this flickering stimulus? 

%%
% To look at this, we calculate the mean average stimulus over some window history length for every time point t, for various different window history lengths. The effect of this operation is shown in the figure below, where the legend refers to the history window in ms. 
history_lengths = [30 102 150 201 300 402 501 600 801 1002];
hl = history_lengths/3; % history lengths better be divisible by 3!
shat = NaN(length(hl),length(PID));
for i = 1:length(hl)
	shat(i,:) = filter(ones(1,hl(i))/hl(i),1,PID);
	shat(i,1:hl(i)) = NaN;
end
figure('outerposition',[10 10 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
plot(time,PID,'k','LineWidth',2)
plot(time,shat(1,:),'r','LineWidth',2)
plot(time,shat(5,:),'b','LineWidth',3)
plot(time,shat(10,:),'g','LineWidth',4)
set(gca,'XLim',[mean(time)-2 mean(time)+2],'LineWidth',2,'box','on','FontSize',font_size)
legend 0ms 30ms 300ms 1002ms 
xlabel('Time (s)','FontSize',font_size)
ylabel('PID (a.u.)','FontSize',font_size)
%%
% Now, we separate neuron responses at times when the mean average stimulus is in the lowest 10% or the highest 10%. These points are marked either green (lowest 10%) or red (or highest 10%) in the figures below, while all the data is plotted in grey. Lines are fit to each of these clouds of points, and the slopes (representing the instantaneous gain) is calculated from these lines. The _y_ axis is the actual firing rate, while the _x_ axis is the predicted firing rate.

% make data vectors
low_slopes = NaN*history_lengths;
low_slopes_err = NaN*history_lengths;
high_slopes = NaN*history_lengths;
high_slopes_err = NaN*history_lengths;
low_gof = NaN*history_lengths;
high_gof = NaN*history_lengths;

% calculate the slopes for all points
[fall gof] = fit(fp(filter_length+2:end)',f(filter_length+2:end)','Poly1');
all_gof = gof.rsquare;
er = confint(fall);
all_slopes_err=fall.p1-er(1,1);
all_slopes = fall.p1;
output_data.all_slopes = all_slopes;

f3= figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
for i = 1:length(history_lengths)

	figure(f3), hold on
	subplot(2,5,i), hold on

	this_shat = shat(i,:);
	this_shat(1:hl(i)) = Inf; % the initial segment where we can't estimate shat is excluded
	[sorted_shat idx] = sort(this_shat,'ascend');
	f_low = f(idx(1:floor(length(PID)/10)));
	fp_low = fp(idx(1:floor(length(PID)/10)));

	this_shat(1:hl(i)) = -Inf;
	[sorted_shat idx] = sort(this_shat,'descend');
	f_high = f(idx(1:floor(length(PID)/10)));
	fp_high = fp(idx(1:floor(length(PID)/10)));

	plot(fp,f,'.','MarkerSize',marker_size,'MarkerFaceColor',[0.7 0.7 0.7],'MarkerEdgeColor',[0.7 0.7 0.7])
	plot(fp_low,f_low,'.','MarkerSize',marker_size,'MarkerFaceColor',[0.5 1 0.5],'MarkerEdgeColor',[0.5 1 0.5])
	plot(fp_high,f_high,'.','MarkerSize',marker_size,'MarkerFaceColor',[1 0.5 0.5],'MarkerEdgeColor',[1 0.5 0.5])
	set(gca,'LineWidth',2,'box','on','FontSize',font_size)
	set(gca,'XLim',[min(f)-10,max(f)+10],'YLim',[min(f)-10,max(f)+10])
	axis square
	title(strcat(mat2str(history_lengths(i)),'ms'))

	if i == 6
		xlabel('Linear Prediction (Hz)','FontSize',16)
		ylabel('Neuron Response (Hz)','FontSize',16)
	end


	% remove NaN values
	f_high(isnan(fp_high)) = [];
	fp_high(isnan(fp_high)) = [];
	f_low(isnan(fp_low)) = [];
	fp_low(isnan(fp_low)) = [];

	% fit lines
	[flow gof] = fit(fp_low',f_low','Poly1');
	low_gof(i) = gof.rsquare;
	er = confint(flow);
	low_slopes_err(i)=flow.p1-er(1,1);
	low_slopes(i) = flow.p1;

	[fhigh gof] = fit(fp_high',f_high','Poly1');
	high_gof(i) =  gof.rsquare;
	er = confint(fhigh);
	high_slopes_err(i)=fhigh.p1-er(1,1);
	high_slopes(i) = fhigh.p1;



	% plot the lines

	plot(min(f):max(f),fall(min(f):max(f)),'Color',[0.5 0.5 0.5],'LineWidth',2)
	plot(min(f_low):max(f_low),flow(min(f_low):max(f_low)),'g','LineWidth',2)
	plot(min(f_high):max(f_high),fhigh(min(f_high):max(f_high)),'r','LineWidth',2)




end



%%
% The following plot shows how the slope of the lines of best fit, or the instantaneous gains, varies with the history length. The plot on the right shows the goodness of fit for each fit, indicating the regions where the fit is meaningful.

% plot to summary figure
figure('outerposition',[10 10 1000 500],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(1,2,1), hold on
plot(history_lengths,all_slopes*ones(1,length(history_lengths)),'k','LineWidth',2), hold on
errorbar(history_lengths,low_slopes,low_slopes_err,'g','LineWidth',2), hold on
errorbar(history_lengths,high_slopes,high_slopes_err,'r','LineWidth',2)
set(gca,'LineWidth',2,'FontSize',20,'box','on','XLim',[0 max(history_lengths)])
xlabel('History Length (ms)','FontSize',20)
ylabel('Slope data/prediction (gain)','FontSize',20)

subplot(1,2,2), hold on
plot(history_lengths,low_gof,'g','LineWidth',2),hold on
plot(history_lengths,high_gof,'r','LineWidth',2)
xlabel('History Length (ms)','FontSize',20)
ylabel('Goodness of Fit, rsquare','FontSize',20)
set(gca,'LineWidth',2,'FontSize',20,'box','on','YLim',[0 1.1],'XLim',[0 max(history_lengths)])

output_data.high_slopes = high_slopes;
output_data.high_slopes_err = high_slopes_err;
output_data.low_slopes = low_slopes;
output_data.low_slopes_err = low_slopes_err;
output_data.low_gof = low_gof;
output_data.high_gof = high_gof;

%% Analysis of Linear Prediction - Filter Variation 
% We have segmented the entire data set based on when the stimulus, averaged in some way over the past, into high and low. Here, I will calculate filters for each of these subsets. 

%% Analysis of Linear Prediction -  Instantaneous Gain 
% There is a mismatch between the linear prediction and the actual firing rate. Moreover, the instantaneous gain seems to be modulated by something that depends on the past history of the stimulus. Here, in the figure below, the instantaneous gain, i.e., the ratio of the actual firing rate to the predicted firing rate, is plotted along with the stimulus. 

gain = f./fp;
figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(2,1,1), hold on
plot(time,gain,'k','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'FontSize',font_size,'LineWidth',2,'box','on','YLim',[0 2])
ylabel('Gain','FontSize',font_size)
subplot(2,1,2), hold on
plot(time,PID,'k','LineWidth',2)
ylabel('PID','FontSize',font_size)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'FontSize',font_size,'LineWidth',2,'box','on')
xlabel('Time (s)','FontSize',font_size)

%%
% In the above analysis, we have considered how a boxcar average over the stimulus immediately preceding the current time affects gain. Now, we want to find some optimal way of averaging the past stimulus history to predict the instantaneous gain: by doing so, we can then predict instantaneous gain, and thus get a better predictor of the actual firing rate.

%%
% In effect, we can calculate a new filter, $K_g$ such that 
% 
% $$ K_g=C\setminus(s'*g) $$
%
% where _g_ is the instantaneous gain and _s_ is the stimulus. 

CRange = [1e-5 1e4];
DRange = [1e-3 1e2];

if ~(exist('KgC') == 1)
	[KgC, Kg,diagnostics_gain] = FindBestFilter(PID(filter_length+2:end),gain(filter_length+2:end),filter_length,CRange,DRange);
end
figure('outerposition',[0 0 400 400],'PaperUnits','points','PaperSize',[900 400]); hold on
plot(filtertime,KgC,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude','FontSize',font_size)
title('Filter: PID > gain','FontSize',20)

%%
% Once again, we can vary the free parameter to see how to best choose a filter. 

PlotFilterDiagnostics(diagnostics_gain);

%%
% From this filter and the stimulus, we can estimate the instantaneous gain (where the instantaneous gain is defined to be the ratio of the actual firing rate to the linear predicted firing rate)
gp = filter(Kg,1,PID - mean(PID)) + 1;
figure('outerposition',[10 10 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
plot(time,gain,'LineWidth',2)
plot(time,gp,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'box','on','LineWidth',2,'FontSize',18)
xlabel('Time','FontSize',20)
ylabel('Gain','FontSize',font_size)
legend Gain PredictedGain

%%
% Using this estimate, we can then correct the linear prediction of the firing rate.
fpg = fp.*gp;
figure('outerposition',[10 10 800 400],'PaperUnits','points','PaperSize',[850 400]); hold on
plot(time,f,'LineWidth',2)
plot(time,fp,'r','LineWidth',2)
plot(time,fpg,'k','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'box','on','LineWidth',2,'FontSize',18)
xlabel('Time','FontSize',20)
ylabel('ORN Firing Rate (Hz)','FontSize',font_size)
legend Data LinearPrediction GainCorrected

%% 
% Is the gain-corrected prediction any better than the linear prediction? The mean error of the gain-corrected prediction is 

disp(sqrt(sum(((fpg(filter_length+2:end)-f(filter_length+2:end)).^2)))*mean(diff(time)))

%%
% while the mean error of the simple linear prediction is
disp(sqrt(sum(((fp(filter_length+2:end)-f(filter_length+2:end)).^2)))*mean(diff(time)))

%%
% No matter what I do, I can't get use a linear filter to predict the gain well from the stimulus. 


%% Summary/Problems
%
% # Best fit line of linear prediction to actual data does not have a slope of 1 for Damon's function. For synthetic data the slope of the line is indeed 1. However, in the real data set, it is not, even with no regularisation. 
% # Prediction of gain from past stimulus is very poor. The correlation coefficient of the prediction with the instantaneous gain is ~ 0.2. For comparison, the corr. coeff. of the firing rate prediction > 0.95. 




%% Next Steps
% 
% # Calculate filters for high and low stimulus zones. 
% # Calculate filters after using elastic net regularisation. 
% # Check if someone has looked at Weber's law in olfaction/ORNs
% # shuffle data to check validity?
% # fit Dynamical Adaptation model to data?

% close all to remove all extraneous figures, since we are publishing to HTML anyway
close all