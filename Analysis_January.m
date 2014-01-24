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

	% make all vectors consistent
	PID = PID(:);
	time = time(:);
	f = f(:);
	fs = fs(:);
	Valve = Valve(:);

	% detrend PID
	ptrend = fit(time,PID,'Poly1'); 
	PID = PID - (ptrend(time) - mean(ptrend(time)));

	% detrend f
	ptrend = fit(time,f,'Poly1'); 
	f = f - (ptrend(time) - mean(ptrend(time)));

	clear ptrend
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
% Details about the filter extraction, regularisation methods, and validation with synthetic and real data is listed in (FilterExtraction.pdf). Here, we calculate the best filter (scaled to ensure that gain = 1) using techniques described in that document from the PID (left) and from the Valve (right).

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
filtertime = 0:mean(diff(time)):filter_length*mean(diff(time)); % this is the x axis for the 

if ~(exist('K') == 1)
	[K diagnostics] = FindBestFilter(PID,f,filter_length);
	K = K*diagnostics.slope(diagnostics.bestfilter);
end

if ~(exist('K_Valve') == 1)
	[K_Valve diagnostics_valve] = FindBestFilter(Valve,f,filter_length);
	K_Valve = K_Valve*diagnostics_valve.slope(diagnostics_valve.bestfilter);
end

% plot
figure('outerposition',[0 0 700 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,K,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude (Hz)','FontSize',font_size)
title('PID > f','FontSize',font_size)

subplot(1,2,2), hold on
plot(filtertime,K_Valve,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude (Hz)','FontSize',font_size)
title('Valve > f','FontSize',font_size)

%%
% The variation of filter shape with regularisation parameter for the PID > f filter is shown below:

PlotFilterDiagnostics2(diagnostics,marker_size,marker_size2,font_size,'PID > f')

%%
% The variation of filter shape with regularisation parameter for the Valve > f filter is shown below:

PlotFilterDiagnostics2(diagnostics_valve,marker_size,marker_size2,font_size,'Valve > f')



%%
% The convolution of this filter with the input (PID) gives a prediction of the output. In the figure below, the data is shown in black, and the prediction from the PID filter is shown in blue and the prediction from the Valve filter is shown in red. Note that sometimes the prediction can be below zero, which doesn't mean anything physical (firing rate cannot be < 0). 

figure('outerposition',[0 0 800 500],'PaperUnits','points','PaperSize',[800 500]); hold on
fp = filter(K,1,PID-mean(PID)) + mean(f);
% censor initial prediction
fp(1:filter_length+1)=NaN;

fp_valve = filter(K_Valve,1,Valve-mean(Valve))+mean(f);
% censor initial prediction
fp_valve(1:filter_length+1)=NaN;

plot(time,f,'k','LineWidth',2)
plot(time,fp,'b','LineWidth',2)
plot(time,fp_valve,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-2 mean(time)+2],'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Time (s)','FontSize',20)
ylabel('Firing Rate (Hz)','FontSize',font_size)
title('Comparison of data and prediction.')
legend Data PID Valve


%% Analysis of Linear Prediction - Which input predicts output better?
% The measurement of the stimulus should help us predict the response of the neuron better than just the valve control signal, as we now capture the fast fluctuations of the odour stimulus as it reaches the ORN. Thus, a prediction of ORN firing response from the PID signal should be better than a prediction of the ORN firing response from the valve signal. 

%%
% The r-square (coefficient of determination) of the PID>f prediction is:
disp(rsquare(f(filter_length+2:end),fp(filter_length+2:end)))

%%
% The r-square (coefficient of determination) of the Valve>f prediction is:
disp(rsquare(f(filter_length+2:end),fp_valve(filter_length+2:end)))


%% Analysis of Linear Prediction - Linearity of Prediction 
% For a linear filter calculated from the data, a plot of the actual response to the predicted response can be fit with a line of slope 1. Here the actual firing rate is plotted on the X axis and the firing rate predicted from the linear filter is plotted on the Y axis. 



figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(fp,f,'.','MarkerSize',marker_size)
set(gca,'LineWidth',2,'FontSize',20,'box','on')
set(gca,'XLim',[min([fp; f]) max([fp; f])],'YLim',[min([fp; f]) max([fp; f])])
axis square
ylabel('Actual Firing rate (Hz)','FontSize',20)
xlabel('Linear Predictionf Firing rate (Hz)','FontSize',20)
[fall ~] = fit(fp(filter_length+2:end),f(filter_length+2:end),'Poly1');
plot(fall(min(f):1:max(f)),min(f):1:max(f),'k','LineWidth',2)
title('PID > f','FontSize',font_size)

subplot(1,2,2), hold on
plot(fp_valve,f,'.','MarkerSize',marker_size)
set(gca,'LineWidth',2,'FontSize',20,'box','on')
set(gca,'XLim',[min([fp_valve; f]) max([fp_valve; f])],'YLim',[min([fp_valve; f]) max([fp_valve; f])])
axis square
xlabel('Linear Prediction Firing rate (Hz)','FontSize',20)
[fall2 ~] = fit(fp_valve(filter_length+2:end),f(filter_length+2:end),'Poly1');
plot(fall2(min(f):1:max(f)),min(f):1:max(f),'k','LineWidth',2)
title('Valve > f','FontSize',font_size)

%% 
% The line is a best fit to all the data. The slope of the best fit line for the prediction from PID is: 
disp(fall.p1)

%% 
% The line is a best fit to all the data. The slope of the best fit line for the prediction from the valve is: 
disp(fall2.p1)


%% Analysis of Linear Prediction - Response to High and Low Stimuli
% How does the response of the ORN differ for high and low stimuli? Specifically, does the neuron display the characteristics of fast adaptation to this flickering stimulus? 

%%
% To look at this, we calculate the mean average stimulus over some window history length for every time point _t_, for various different window history lengths. The effect of this operation is shown in the figure below, where the legend refers to the history window in ms. 
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
% Now, we separate neuron responses and linear prediction at times when the mean average stimulus is in the lowest 10% or the highest 10%. These points are marked either green (lowest 10%) or red (or highest 10%) in the figures below, while all the data is plotted in grey. Lines are fit to each of these clouds of points, and the slopes (representing the instantaneous gain) is calculated from these lines. The _y_ axis is the actual firing rate, while the _x_ axis is the predicted firing rate. In this analysis the prediction from the PID is used.

% make data vectors
low_slopes = NaN*history_lengths;
low_slopes_err = NaN*history_lengths;
high_slopes = NaN*history_lengths;
high_slopes_err = NaN*history_lengths;
low_gof = NaN*history_lengths;
high_gof = NaN*history_lengths;

% calculate the slopes for all points
[fall gof] = fit(fp(filter_length+2:end),f(filter_length+2:end),'Poly1');
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
	[flow gof] = fit(fp_low,f_low,'Poly1');
	low_gof(i) = gof.rsquare;
	er = confint(flow);
	low_slopes_err(i)=flow.p1-er(1,1);
	low_slopes(i) = flow.p1;

	[fhigh gof] = fit(fp_high,f_high,'Poly1');
	high_gof(i) =  gof.rsquare;
	er = confint(fhigh);
	high_slopes_err(i)=fhigh.p1-er(1,1);
	high_slopes(i) = fhigh.p1;


	% plot the best fit lines
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

%% Analysis of Linear Prediction - Filter Variation (to be done)
% We have segmented the entire data set based on when the stimulus, averaged in some way over the past, into high and low. Here, I will calculate filters for each of these subsets. 

%% Analysis of Linear Prediction - Does adding another linear filter improve prediction?
% An ideal linear filter should capture all the linear variation in the data. This means that if one were to construct a vector of residuals, a linear filter would be unable to predict the residuals from the data, and would lead to no improvement on the original filter. Is this true? 

%%
% First, we construct a vector of residuals by subtracting the linear prediction from the data.
res = f - fp;
figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(2,1,1), hold on
plot(time,res,'k','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'FontSize',font_size,'LineWidth',2,'box','on')
ylabel('Residuals','FontSize',font_size)
subplot(2,1,2), hold on
plot(time,PID,'k','LineWidth',2)
ylabel('PID','FontSize',font_size)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'FontSize',font_size,'LineWidth',2,'box','on')
xlabel('Time (s)','FontSize',font_size)

%%
% we then construct a filter from the PID to these residuals to try to predict the residuals.
[Kres ,diagnostics_res] = FindBestFilter(PID(filter_length+2:end),res(filter_length+2:end),filter_length);

figure('outerposition',[0 0 400 400],'PaperUnits','points','PaperSize',[800 400]); hold on
plot(filtertime,Kres,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude','FontSize',font_size)
title('Filter: PID > residuals','FontSize',20)

%%
% and use this filter to predict the residuals from the stimulus. 
resp = filter(Kres,1,PID-mean(PID)) + mean2(res);
figure('outerposition',[10 10 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
plot(time,res,'k','LineWidth',2)
plot(time,resp,'b','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'box','on','LineWidth',2,'FontSize',18)
xlabel('Time','FontSize',20)
ylabel('Residuals','FontSize',font_size)
legend Residuals 'Predicted Residuals'

%%
% Does adding this back to the prediction improve the prediction? The following plot shows the data (black), the linear prediction (red), and the linear prediction corrected by the prediciton of the residual (green).
fpr = fp + resp;
figure('outerposition',[10 10 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
plot(time,fpr,'g','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'box','on','LineWidth',2,'FontSize',18)
xlabel('Time','FontSize',20)
ylabel('Residuals','FontSize',font_size)
legend Data 'Linear Prediction' 'Residual Corrected'

%%
% The r-square of the corrected prediction is
disp(rsquare(f(filter_length+2:end),fpr(filter_length+2:end)))

%% 
% and the r-square of the simple linear prediction is
disp(rsquare(f(filter_length+2:end),fp(filter_length+2:end)))

%%
% How does this make any sense? If this is true, then simply adding the filters together will yield a better filter, which we should have computed in the very beginning. The filter, the second filter computed from residuals, and their sum is shown in the panel on the left. On the right is a filter with lower regularisation. 

figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,K,'k','LineWidth',2)
plot(filtertime,Kres,'r','LineWidth',2)
plot(filtertime,Kres+K,'g','LineWidth',2)
legend Filter ResidualFilter Sum
set(gca,'FontSize',font_size,'LineWidth',2,'box','on')
ylabel('Filter Amplitude (Hz)','FontSize',font_size)

K_lowreg = FitFilter2Data(PID,f,'reg=1e-1;');
subplot(1,2,2), hold on
plot(filtertime,K_lowreg,'k','LineWidth',2)
title('Filter with low regularisation','FontSize',font_size)
set(gca,'FontSize',font_size,'LineWidth',2,'box','on')
xlabel('Time (s)','FontSize',font_size)

%%
% So what we are doing is simply undoing the work we did in regularising the filter. 

%% Analysis of Gain: Instantaneous Gain and Fitted Gain
% There is a mismatch between the linear prediction and the actual firing rate. Moreover, the instantaneous gain seems to be modulated by something that depends on the past history of the stimulus. Here, in the figure below, the instantaneous gain, i.e., the ratio of the actual firing rate to the predicted firing rate, is plotted along with the stimulus. 

g = f./fp;
figure('outerposition',[0 0 800 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(2,1,1), hold on
plot(time,g,'k','LineWidth',2)
set(gca,'XLim',[mean(time)-2 mean(time)+2],'FontSize',font_size,'LineWidth',2,'box','on','YLim',[0 2])
ylabel('Gain','FontSize',font_size)
subplot(2,1,2), hold on
plot(time,PID,'k','LineWidth',2)
ylabel('PID','FontSize',font_size)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'FontSize',font_size,'LineWidth',2,'box','on')
xlabel('Time (s)','FontSize',font_size)

%%
% We can also calculate the moving gain by fitting lines to the data and the prediction collected in bins of size _w_. The following plots show the gain computed in this way for a few different _w_. We want to do this because the raw instantaneous gain is very noisy, and predicting it from PID is hard. 

% WARNING. THIS CODE IS SUPER SLOW. FITTING THOUSANDS OF LINES. 
ws = [5 10 30 100];
if ~(exist('gain') == 1)
	gain = f./fp;
	gain = repmat(gain,1,length(ws)+1);
	for i = 1:length(ws)
		disp(i)
		starthere = filter_length+2+ws(i);
		for j = starthere:1:length(f)
			thesepoints = (j - ws(i)):j;
			fitmetrics = fit(fp(thesepoints),f(thesepoints),'Poly1');
			gain(j,i+1) = fitmetrics.p1;
		end
	end
end

ws = [1 ws];
figure('outerposition',[10 10 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold all
legendtext = {};
for i = 1:size(gain,2)
	plot(time,gain(:,i))
	legendtext{i} = strcat('Window size:',mat2str(ws(i)*3),'ms');
end
set(gca,'XLim',[mean(time)-2 mean(time)+2],'LineWidth',2,'box','on','FontSize',font_size)
legend(legendtext)
xlabel('Time (s)','FontSize',font_size)
ylabel('Gain (f/fp)','FontSize',font_size)

%%
% Of these different gain vectors, which fixes the prediction the best? The following vector is the r-square of prediction, corrected by the each of the gain vectors, scaled by the rsquare of the simple linear prediciton. The first number is > 1, indicating that it is an improvement over the simple linear prediction. (In fact, by definition, the first number corresponds to a perfect prediction)
fpg = f;
fpg = repmat(fpg,1,length(ws));
for i = 1:length(ws)
	fpg(:,i) = fp.*gain(:,i);
end

rvalues = NaN(1,length(ws));
for i = 1:length(ws)
	rvalues(i) = rsquare(f(filter_length+2:end),fpg(filter_length+2:end,i));
end
disp(rvalues/(rsquare(f(filter_length+2:end),fp(filter_length+2:end))));

%%
% Note that none of the scaled r-square values are > 1 for any of the fits (points > 3ms) so this method of finding the gain is not helpful in improving the prediction. 

%% Analysis of Gain: Smoothed gain. 
% Instead of calculating the instantaneous point-by-point gain, we can also smooth the instantaneous gain over some small window size _w_. The following plot shows the effect of smoothing for a few window sizes _w_. 
gain2 = f./fp;
ws = [1 5 10 25 50 100];
gain2 = repmat(gain2,1,length(ws));
for i = 1:length(ws)
	gain2(:,i)=smooth(gain2(:,i),ws(i));	
end

figure('outerposition',[10 10 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold all
legendtext = {};
for i = 1:size(gain2,2)
	plot(time,gain2(:,i),'LineWidth',2)
	legendtext{i} = strcat('Window size:',mat2str(ws(i)*3),'ms');
end
set(gca,'XLim',[mean(time)-2 mean(time)+2],'LineWidth',2,'box','on','FontSize',font_size)
legend(legendtext)
xlabel('Time (s)','FontSize',font_size)
ylabel('Gain (f/fp)','FontSize',font_size)

%%
% How does smoothing the gain affect our ability to correct the linear prediction? The following plot shows how the r-square of the corrected linear prediction varies with smoothing the gain. Also shown is the line which indicates the linear prediction of the uncorrected simple linear prediction. 

fpg = f;
fpg = repmat(fpg,1,length(ws));
for i = 1:length(ws)
	fpg(:,i) = fp.*gain2(:,i);
end

rvalues = NaN(1,length(ws));
for i = 1:length(ws)
	rvalues(i) = rsquare(f(filter_length+2:end),fpg(filter_length+2:end,i));
end

figure('outerposition',[10 10 400 400],'PaperUnits','points','PaperSize',[400 400]); hold on
plot(ws*3,rvalues,'.')
xlabel('Smoothing window (ms)','FontSize',font_size)
ylabel('r-square of corrected fit','FontSize',font_size)
line([1 max(ws)*3],[rsquare(f(filter_length+2:end),fp(filter_length+2:end)) rsquare(f(filter_length+2:end),fp(filter_length+2:end))],'LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size)



%%
% In the above analysis, we have considered how a boxcar average over the stimulus immediately preceding the current time affects gain. Now, we want to find some optimal way of averaging the past stimulus history to predict the instantaneous gain: by doing so, we can then predict instantaneous gain, and thus get a better predictor of the actual firing rate.

%%
% In effect, we can calculate a new filter, $K_g$ such that 
% 
% $$ K_g=\hat{C}\setminus(s'*g) $$
%
% where _g_ is the instantaneous gain (the point-by-point ratio of the data to the linear prediction) and _s_ is the stimulus. 


if ~(exist('Kg_PID') == 1)
	[Kg_PID ,diagnostics_gainPID] = FindBestFilter(PID(filter_length+2:end),gain(filter_length+2:end),filter_length);
end

if ~(exist('Kg_Valve') == 1)
	[Kg_Valve ,diagnostics_gainv] = FindBestFilter(f(filter_length+2:end),gain(filter_length+2:end),filter_length);
end

figure('outerposition',[0 0 800 400],'PaperUnits','points','PaperSize',[800 400]); hold on
subplot(1,2,1), hold on
plot(filtertime,Kg_PID,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude','FontSize',font_size)
title('Filter: PID > gain','FontSize',20)

subplot(1,2,2), hold on
plot(filtertime,Kg_Valve,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size)
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude','FontSize',font_size)
title('Filter: Valve > gain','FontSize',20)

%%
% Once again, we can look at how filter shape is affected by our choise of _r_ for the PID > gain prediction: 
PlotFilterDiagnostics2(diagnostics_gainPID,marker_size,marker_size2,font_size,'PID > Gain');

%%
% Once again, we can look at how filter shape is affected by our choise of _r_ for the PID > gain prediction: 
PlotFilterDiagnostics2(diagnostics_gainv,marker_size,marker_size2,font_size,'Valve > Gain');


%%
% From this filter and the stimulus, we can estimate the instantaneous gain (where the instantaneous gain is defined to be the ratio of the actual firing rate to the linear predicted firing rate)
gp_PID = filter(Kg_PID,1,PID - mean(PID)) + 1;
gp_v = filter(Kg_Valve,1,f - mean(f)) + 1;
figure('outerposition',[10 10 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
plot(time,gain,'k','LineWidth',2)
plot(time,gp_PID,'b','LineWidth',2)
plot(time,gp_v,'r','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'box','on','LineWidth',2,'FontSize',18)
xlabel('Time','FontSize',20)
ylabel('Gain','FontSize',font_size)
legend Gain 'PID Prediction' 'Valve Prediction'

%%
% Using these estimates, we can then correct the linear prediction of the firing rate.
fp_gpPID = fp.*gp_PID;
fp_gpv = fp.*gp_v;
figure('outerposition',[10 10 800 400],'PaperUnits','points','PaperSize',[850 400]); hold on
plot(time,f,'k','LineWidth',2)
plot(time,fp,'r','LineWidth',2)
plot(time,fp_gpPID,'b','LineWidth',2)
plot(time,fp_gpv,'g','LineWidth',2)
set(gca,'XLim',[mean(time)-1 mean(time)+2],'box','on','LineWidth',2,'FontSize',18)
xlabel('Time','FontSize',20)
ylabel('ORN Firing Rate (Hz)','FontSize',font_size)
legend Data LinearPrediction 'Gain Corrected (PID)' 'Gain Corrected (Valve)'

%% 
% Is the gain-corrected prediction any better than the linear prediction? The r-square of the PID gain corrected estimate is: 
disp(rsquare(f(filter_length+2:end),fp_gpPID(filter_length+2:end)))

%% 
% The r-square of the Valve gain corrected estimate is: 
disp(rsquare(f(filter_length+2:end),fp_gpv(filter_length+2:end)))

%% 
% This is in comparison to the simple linear prediction with a r-square of:
disp(rsquare(f(filter_length+2:end),fp(filter_length+2:end)))


%%
% The autocorrelation function of the gain shows that it is very tightly constrained, almost as much as the valve, and much less than the PID, which we are trying to use to predict it. 
cf = xcorr(f(filter_length+2:end)-mean(f),'unbiased'); cf=cf/max(cf((length(cf)+1)/2-100:(length(cf)+1)/2+100));
cp = xcorr(PID(filter_length+2:end)-mean(PID),'unbiased'); cp=cp/max(cp((length(cp)+1)/2-100:(length(cp)+1)/2+100));
cv = xcorr(Valve(filter_length+2:end)-mean(Valve),'unbiased'); cv=cv/max(cv((length(cv)+1)/2-100:(length(cv)+1)/2+100));
cg = xcorr(gain(filter_length+2:end)-mean(gain(filter_length+2:end)),'unbiased');  cg=cg/max(cg((length(cg)+1)/2-100:(length(cg)+1)/2+100));
cn2 = xcorr(randn(1,length(f(filter_length+2:end))),'unbiased'); cn2 = cn2/max(cn2);
ctime =  mean(diff(time)):mean(diff(time)):length(cf)*mean(diff(time)); 
ctime = ctime - mean(ctime);
figure('outerposition',[10 10 850 400],'PaperUnits','points','PaperSize',[850 400]); hold on
subplot(1,2,1), hold on
plot(ctime,cf,'b','LineWidth',2)
plot(ctime,cp,'r','LineWidth',2)
plot(ctime,cv,'g','LineWidth',2)
plot(ctime,cg,'k','LineWidth',2)
set(gca,'box','on','XLim',[-0.5 0.5],'FontSize',font_size,'LineWidth',2,'YLim',[0 1.5])
xlabel('Time (s)','FontSize',font_size)
ylabel('Normalised autocorrelation','FontSize',font_size)
legend ORN PID Valve Gain 
subplot(1,2,2), hold on
plot(ctime,cn2,'r','LineWidth',2)
set(gca,'box','on','XLim',[-0.5 0.5],'FontSize',font_size,'LineWidth',2,'YLim',[0 1.5])
xlabel('Time (s)','FontSize',font_size)
legend GaussianNoise



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