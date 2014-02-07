%% Dynamical Adaptation in ORNs 
% Do ORNs exhibit fast adaptation to a flickering stimulus? Can a simple dynamical adaptation model predict ORN responses to flickering inputs? Here, I take data from Carlotta's random stimulus experiments and first check how well a simple linear prediction from the stimulus compares to the data. Then, I study how the instantaneous gain of the actual ORN response w.r.t to the predicted ORN response varies with time, and try to find a objective filter from the stimulus to this instantaneous gain to correct for systematic errors in the linear prediction.  

%%

% parameters to tune for figure display
font_size = 20;
font_size2 = 18;
marker_size = 10;
marker_size2 = 20;

%% Data Visualisation 
% Data from this file will be used for this analysis.

filename = '/data/random-stim/final_2011_06_14_ab3A_1o3ol3X-3_20ml_30sec_30ms_rand.mat';
disp(filename)

%%
% Here, data from simultaneous measurements of the stimulus and from the ORN is shown. The top row shows the valve command signal (a binary signal, high means valve is open). The middle row shows the simultaneous measurement of the odour stimulus with a PID, and the bottom row shows the instantaneous firing rate of the ORN recorded from. In blue is the data computed using Carlotta's code from the raw data. In red is the data computed using my code. Notably different is the firing rate curve, which is much smoother and loses all the high-frequency wiggles. This is because I use a Gaussian smoothing window to do this, while Carlotta uses a boxcar window. 

crop_traces = 0;
% grab data
if ~(exist('PID') == 1)
	crop_traces = 1;
	[PID1,time1, f1, ~, Valve] = PrepData2(filename);
	[PID2, time2, f2, Valve2] = PrepData3(filename);

	% make all vectors consistent
	PID1 = PID1(:);
	time1 = time1(:);
	f1 = f1(:);
	Valve = Valve(:);
	PID2 = PID2(:);
	time2 = time2(:);
	f2 = f2(:);
	Valve2 = Valve2(:);

	% correct offset
	PID2 = PID2-min(PID2);
	PID2 = PID2/max(PID2);
	PID2 = PID2+min(PID1);
	PID2 = PID2/max(PID2);
	PID2 = PID2*max(PID1);




	% detrend PID
	ptrend = fit(time1,PID1,'Poly1'); 
	PID1 = PID1 - (ptrend(time1) - mean(ptrend(time1)));

	% detrend f
	ptrend = fit(time1,f1,'Poly1'); 
	f1 = f1 - (ptrend(time1) - mean(ptrend(time1)));

	ptrend = fit(time2,PID2,'Poly1'); 
	PID2 = PID2 - (ptrend(time2) - mean(ptrend(time2)));

	% detrend f
	ptrend = fit(time2,f2,'Poly1'); 
	f2 = f2 - (ptrend(time2) - mean(ptrend(time2)));

	clear ptrend
end



% plot the data
figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(3,1,1), hold on
plot(time1,Valve,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'XLim',[18 20],'YLim',[-0.1 1.1])
ylabel('stimulus')

subplot(3,1,2), hold on
plot(time1,PID1,'b','LineWidth',1.5)
plot(time2,PID2,'r','LineWidth',1.5)
set(gca,'box','on','LineWidth',2,'XLim',[18 20])
ylabel('PID (a.u.)')

subplot(3,1,3), hold on
plot(time1,f1,'b','LineWidth',1.5)
plot(time2,f2,'r','LineWidth',1.5)
set(gca,'box','on','LineWidth',2,'XLim',[18 20])
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

% use my estimates
PID = PID2;
f = f2;
time = time1;



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
end


% compute filter
filter_length = 333;
filtertime = 0:mean(diff(time)):filter_length*mean(diff(time));  % this is the time axis for the filter.
filtertime = filtertime - shift_input*(mean(diff(time))); % correct for shifted input

if ~(exist('K') == 1)
	[K diagnostics] = FindBestFilter(PID,f);
end

if ~(exist('K_Valve') == 1)
	[K_Valve diagnostics_valve] = FindBestFilter(Valve,f);
end

% plot
figure('outerposition',[0 0 700 350],'PaperUnits','points','PaperSize',[800 350]); hold on
subplot(1,2,1), hold on
plot(filtertime,K,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size,'XLim',[min(filtertime) max(filtertime)])
xlabel('Lag (s)','FontSize',20)
ylabel('Filter Amplitude (Hz)','FontSize',font_size)
title('PID > f','FontSize',font_size)

subplot(1,2,2), hold on
plot(filtertime,K_Valve,'LineWidth',2)
set(gca,'box','on','LineWidth',2,'FontSize',font_size,'XLim',[min(filtertime) max(filtertime)])
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
% The lines are best fits to all the data. The slopes of these lines should be exactly 1 as the "best filters" are supposed to be unit gain. The slope, in fact is:
disp(fall.p1)



%% Analysis of Linear Prediction - Response to High and Low Stimuli
% How does the response of the ORN differ for high and low stimuli? Specifically, does the neuron display the characteristics of fast adaptation to this flickering stimulus? 

%%
% To look at this, we calculate the mean average stimulus over some window history length for every time point _t_, for various different window history lengths. The effect of this operation is shown in the figure below, where the legend refers to the history window in ms. 

%%
% To calcualte the mean average stimulus, I use the _filtfilt_ function, which is a zero-phase digital filter that works by passing it through once in the forward direction, and once again in the reverse direction. I previously used the simple _filter_ function, but the results are substantially different: now, the separation of gains extends for larger history window lengths. 

history_lengths = [30 102 150 201 300 402 501 600 801 1002];
hl = history_lengths/3; % history lengths better be divisible by 3!
shat = NaN(length(hl),length(PID));
for i = 1:length(hl)
	shat(i,:) = filtfilt(ones(1,hl(i))/hl(i),1,PID);
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

[output_data] = GainAnalysis(f,fp,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,1);
%%
% The following plot shows how the slope of the lines of best fit, or the instantaneous gains, varies with the history length. The plot on the right shows the goodness of fit for each fit, indicating the regions where the fit is meaningful.

[output_data] = GainAnalysis(f,fp,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,[2 3]);

%% Analysis of Linear Prediction - Filter Variation 
% We have segmented the entire data set based on when the stimulus, averaged in some way over the past, into high and low. The following plots show how the filter shape for low, high and middle 10% of the stimuli changes with the history window size. In these plots, only 10% of the data, either the lowest 10% (green), the highest 10% (red) of the middle 10% (black), based on the window smoothing is used to compute each of the filters.


Kmiddle = NaN(length(history_lengths),filter_length+1);
Klow = NaN(length(history_lengths),filter_length+1);
Khigh = NaN(length(history_lengths),filter_length+1);
for i = 1:length(history_lengths)

	slicelength = floor(length(PID)/10);
	% low filter
	this_shat = shat(i,:);
	this_shat(1:hl(i)) = Inf; % the initial segment where we can't estimate shat is excluded
	[sorted_shat idx] = sort(this_shat,'ascend');
	OnlyThesePoints=idx(1:slicelength);
	Klow(i,:)=FitFilter2Data(PID,f,OnlyThesePoints);

	% low filter
	this_shat = shat(i,:);
	this_shat(1:hl(i)) = Inf; % the initial segment where we can't estimate shat is excluded
	[sorted_shat idx] = sort(this_shat,'descend');
	OnlyThesePoints=idx(1:slicelength);
	Khigh(i,:)=FitFilter2Data(PID,f,OnlyThesePoints);

	% middle filter
	OnlyThesePoints=idx((length(idx)/2-slicelength/2):(length(idx)/2+slicelength/2));
	Kmiddle(i,:)=FitFilter2Data(PID,f,OnlyThesePoints);

	
end

% plot
f3= figure('outerposition',[0 0 1200 600],'PaperUnits','points','PaperSize',[1200 600]); hold on
for i = 1:length(history_lengths)
	subplot(2,5,i), hold on

	plot(filtertime,Kmiddle(i,:),'k','LineWidth',2)
	plot(filtertime,Klow(i,:),'g','LineWidth',2)
	plot(filtertime,Khigh(i,:),'r','LineWidth',2)
	set(gca,'FontSize',font_size2,'LineWidth',2,'YLim',[min(min(Klow))-5 max(max(Klow))+5],'box','on','XLim',[min(filtertime) max(filtertime)])
	title(strcat(mat2str(history_lengths(i)),'ms'))

	if i == 6
		xlabel('Filter Lag (s) ','FontSize',16)
		ylabel('Neuron Response (Hz)','FontSize',16)
	end
end


%%
% We can also plot the height of the filter _vs._ the history window size:

figure('outerposition',[10 10 800 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(history_lengths,max(Kmiddle'),'k','LineWidth',2)
plot(history_lengths,max(Klow'),'g','LineWidth',2)
plot(history_lengths,max(Khigh'),'r','LineWidth',2)
set(gca,'FontSize',font_size,'LineWidth',2,'YLim',[0 max(max(Klow))+5],'box','on')
xlabel('History Length (ms) ','FontSize',font_size)
ylabel('Filter Height (Hz)','FontSize',font_size)

%%
% We also expect the low stimulus (high gain) filter to be slower than the high stimulus (low gain) filter. Is this so in the data? Here, I plot the time of peak of each filter _vs._ the history window size. 
figure('outerposition',[10 10 800 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
[~,tmax]=max(Kmiddle'); tmax = tmax*mean(diff(time)) - shift_input*mean(diff(time));
plot(history_lengths,tmax,'k','LineWidth',2)
[~,tmax]=max(Klow'); tmax = tmax*mean(diff(time)) - shift_input*mean(diff(time));
plot(history_lengths,tmax,'g','LineWidth',2)
[~,tmax]=max(Khigh'); tmax = tmax*mean(diff(time)) - shift_input*mean(diff(time));
plot(history_lengths,tmax,'r','LineWidth',2)
set(gca,'FontSize',font_size,'LineWidth',2,'YLim',[0 2*max(max(tmax))],'box','on')
xlabel('History Length (ms) ','FontSize',font_size)
ylabel('Time of filter peak (s)','FontSize',font_size)


%% Analysis of Gain: Instantaneous Gain and Fitted Gain
% There is a mismatch between the linear prediction and the actual firing rate. Moreover, the instantaneous gain seems to be modulated by something that depends on the past history of the stimulus (see analysis in the previous section). We also know that we cannot predict with a simple linear filter any linear errors that the filter makes. However, there is the possibility of fitting a linear filter to the gain, which is a non-linear function of the filter output. But how do we compute this gain? Here, in the figure below, the instantaneous gain, i.e., the ratio of the actual firing rate to the predicted firing rate, is plotted along with the stimulus. 

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
% we're going to cheat and pre-load data
ws = [5 10 30 100];
if ~(exist('gain') == 1)
	% gain = f./fp;
	% gain = repmat(gain,1,length(ws)+1);
	% for i = 1:length(ws)
	% 	disp(i)
	% 	starthere = filter_length+2+ws(i);
	% 	for j = starthere:1:length(f)
	% 		j
	% 		fitmetrics = fit(fp(thesepoints),f(thesepoints),'Poly1');
	% 		gain(j,i+1) = fitmetrics.p1;
	% 	end
	% end
	load('gain.mat')
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
% Of these different gain vectors, which fixes the prediction the best? The following plot shows the r-square values of corrected linear prediction, corrected by the gain computed in different ways (lines fit to windows of different lengths). On the y-axis is the r-square, and the x-axis is the window over which the lines are fit.  The first point has a r-square of 1, indicating a perfect fit (this is by definition, since the first point is the instantaneous gain). Note that none of the other points exceed the line, which indicates the r-square of the simple linear prediction. 
fpg = f;
fpg = repmat(fpg,1,length(ws));
for i = 1:length(ws)
	fpg(:,i) = fp.*gain(:,i);
end

rvalues = NaN(1,length(ws));
for i = 1:length(ws)
	rvalues(i) = rsquare(f(filter_length+2:end),fpg(filter_length+2:end,i));
end

figure('outerposition',[10 10 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(ws*3,rvalues,'.','MarkerSize',marker_size2)
xlabel('Window over which slope is fitted (ms)','FontSize',font_size)
ylabel('r-square of corrected fit','FontSize',font_size)
line([1 max(ws)*3],[rsquare(f(filter_length+2:end),fp(filter_length+2:end)) rsquare(f(filter_length+2:end),fp(filter_length+2:end))],'LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size)

%%
% We conclude that this method of finding the gain is not helpful in improving the prediction. 

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

figure('outerposition',[10 10 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(ws*3,rvalues,'.','MarkerSize',marker_size2)
xlabel('Smoothing window (ms)','FontSize',font_size)
ylabel('r-square of corrected fit','FontSize',font_size)
line([1 max(ws)*3],[rsquare(f(filter_length+2:end),fp(filter_length+2:end)) rsquare(f(filter_length+2:end),fp(filter_length+2:end))],'LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size)

%%
% From this analysis of the gain, it is clear that if we do want to improve the linear prediction, the best way to do it (apart from the instantaneous gain, which would lead to a perfect prediction), is smoothing the gain in some way. In the next section, we will try to estimate the smoothed gain.

%% Gain Analysis - Estimating the gain 
% In the above analysis, we have considered how a boxcar average over the stimulus immediately preceding the current time affects gain. Now, we want to find some optimal way of averaging the past stimulus history to predict the instantaneous gain: by doing so, we can then predict gain, and thus get a better predictor of the actual firing rate.

%%
% In effect, we can calculate a new filter, $K_g$ such that 
% 
% $$ K_g=\hat{C}\setminus(s'*g) $$
%
% where _g_ is the smoothed gain and _s_ is the stimulus. 

%%
% For with various smoothing of the instantaneous gain, we find the best filter for each and plot it below:

if ~(exist('Kg_smooth') == 1)
	Kg_smooth = zeros(filter_length+1,size(gain2,2));
	for i = 1:size(gain2,2)
		Kg_smooth(:,i) = FindBestFilter(PID(filter_length+2:end),gain2(filter_length+2:end,i));
	end

end


figure('outerposition',[10 10 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
plot(filtertime,Kg_smooth,'LineWidth',2)
xlabel('Filter Lag','FontSize',font_size)
ylabel('Filter Amplitude','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[min(filtertime) max(filtertime)])
legend(legendtext)

%%
% How well can we predict the smoothed gain vectors using these filters and the stimulus? The following plot shows the r-square between prediction of the smoothed gain and the smoothed gain for various smoothing windows. 

gain2p = gain2;
for i = 1:length(ws)
	gain2p(:,i) = filter(Kg_smooth(:,i),1,PID-mean(PID)) + mean2(gain2(filter_length+2:end,i));
end

rvalues = NaN(1,length(ws));
for i = 1:length(ws)
	rvalues(i) = rsquare(gain2p(filter_length+2:end,i),gain2(filter_length+2:end,i));
end

figure('outerposition',[10 10 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(ws*3,rvalues,'.','MarkerSize',marker_size2)
xlabel('Smoothing window (ms)','FontSize',font_size)
ylabel('r-square','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size)

%%
% Well, that sucks. Maybe even this horrible estimation of the gain can improve the linear prediction of the firing rate in some way? Here, I plot the r-square of the gain-corrected linear prediction (corrected by the prediction of the gain), as a function of gain smoothing (in blue). Also plotted is the r-square of the gain-corrected linear prediction, corrected now with actual gain (in black).


fpg = f;
fpg = repmat(fpg,1,length(ws));
for i = 1:length(ws)
	fpg(:,i) = fp.*gain2(:,i);
end

rvalues = NaN(1,length(ws));
for i = 1:length(ws)
	rvalues(i) = rsquare(f(filter_length+2:end),fpg(filter_length+2:end,i));
end


fpg = f;
fpg = repmat(fpg,1,length(ws));
for i = 1:length(ws)
	fpg(:,i) = fp.*gain2p(:,i);
end

rvalues2 = NaN(1,length(ws));
for i = 1:length(ws)
	rvalues2(i) = rsquare(f(filter_length+2:end),fpg(filter_length+2:end,i));
end

figure('outerposition',[10 10 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(ws*3,rvalues,'b.','MarkerSize',marker_size2), hold on
plot(ws*3,rvalues2,'k.','MarkerSize',marker_size2)
xlabel('Smoothing window (ms)','FontSize',font_size)
ylabel('r-square of corrected fit','FontSize',font_size)
line([1 max(ws)*3],[rsquare(f(filter_length+2:end),fp(filter_length+2:end)) rsquare(f(filter_length+2:end),fp(filter_length+2:end))],'LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size)



%% Re-evaluating filter performance after filtering ORN responses
% The ORN firing rate traces are quite noisy, and it is possible that this is screwing up our estimation of the gain. It is also possible that the linear filter is doing a really good job and can't really be improved, and all the gain residuals are just noise. To check for these possibilities, I will filter the ORN response with a small sliding window, and then recalculate everything. 
fsmooth = filtfilt(ones(1,10)/10,1,f);

% find best filter from PID > smooth f
Kfsmooth = FindBestFilter(PID,fsmooth);
fpsmooth = filter(Kfsmooth,1,PID-mean(PID))+mean(f);
fpsmooth(1:filter_length) = NaN;

figure('outerposition',[10 10 800 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(time,f,'k','LineWidth',1.5)
plot(time,fsmooth,'r','LineWidth',1.5)
plot(time,fpsmooth,'g','LineWidth',1.5)
set(gca,'XLim',[20 22],'FontSize',font_size,'LineWidth',2)
legend Data SmoothedData Prediction

%%
% The r-square of the smoothed ORN response to the actual data is
disp(rsquare(f,fsmooth))

%%
% and the r-square of the smoothed ORN response to the predicted smoothed data is
disp(rsquare(fpsmooth,fsmooth))

%%
% The smoothing operation should not interfere with the gain modulation we observe. Here, we repeat the gain analysis, but this time with the smoothed data and the smoothed prediction. 

[output_data] = GainAnalysis(fsmooth,fpsmooth,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,2);

%%
% Once again, we can estimate the instantaneous gain and construct a filter from the PID to the instantaneous gain. This filter looks like

figure('outerposition',[10 10 400 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
gs = fpsmooth./fsmooth;
Kgsmooth = FindBestFilter(PID(filter_length+2:end),gs(filter_length+2:end));
plot(filtertime,Kgsmooth,'k','LineWidth',2)
title('PID > Instantaneous gain','FontSize',font_size)
xlabel('Filter Lag (s)','FontSize',font_size)
ylabel('Filter Amplitude','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[min(filtertime) max(filtertime)])

% now predict
gsp = filter(Kgsmooth,1,PID-mean(PID)) + 1;
fpg = gsp.*fp;

%%
% This filter can be used to correct the simple linear prediction from the original data. Now, we can compare the gain-corrected prediction to the _original unfiltered firing rate_ and see if is better than the simple linear prediction. 

figure('outerposition',[10 10 800 400],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(time,f,'k','LineWidth',1.5)
plot(time,fp,'b','LineWidth',1.5)
plot(time,fpg,'r','LineWidth',1.5)
set(gca,'XLim',[20 22],'FontSize',font_size,'LineWidth',2)
legend({'Data','Linear Prediction','Gain Corrected'})
xlabel('Time (s)','FontSize',font_size)
ylabel('Firing Rate (Hz)','FontSize',font_size)

%%
% The r-square of the gain-corrected prediction to the data is 
disp(rsquare(fpg,f))

%% 
% _cf._ the  r-square of the simple linear prediction to the data which is 
disp(rsquare(fp,f))

%% 
% Does this gain correction make the gains to low and high stimuli the same? In the following plot, we plot the best-fit lines of corrected prediction to data for different history lengths, and also show an example scatter plot for the gain corrected prediction. 
[output_data] = GainAnalysis(f,fpg,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,2);

h=subplot(1,2,2);


[output_data] = GainAnalysis(f,fpg,PID,shat(history_lengths == 501,:),501,501/3,filter_length,marker_size,marker_size2,font_size,1,h);
xlabel('ORN firing rate (Hz)','FontSize',font_size)
ylabel('Gain corrected prediction (Hz)','FontSize',font_size)


%% A note on elastic net regularisation 
% Carlotta says in her paper that she uses elastic net regularisation (using some unknown, third party code) to minimise the L-1 and L-2 norms of the covariance matrix. However, her code shows that she uses simply L-2 minimisation as discussed previously. I talked to her, and she said that she arbitrarily chose parameters for regularisation anyway, and it is possible that she was using simple L-2 regularisation, as this is a subcase of the elastic net regularisation. 

%% Next Steps
% 
% # Check if Weber's Law is being followed using Carlotta's step responses data
% # Run analysis on all data files. 
% # Fit DA model to data with tau set to zero


%% Docs
% This document was generated by MATLAB's _publish_ function. All files needed to generate this document are on a git repository. To clone onto your machine, use:
%
%   git clone https://srinivasgs@bitbucket.com/srinivasgs/da.git
% 
% You will also need a bunch of functions that this depends on, which are also on git. Use
% 
%   git clone https://srinivasgs@bitbucket.com/srinivasgs/core.git
% 
% to get these. 
% 
% Once you have everything, run these commands to generate this document: 
% 
%   options = struct('showCode',false,'format','latex','imageFormat',...
%   'pdf','figureSnapMethod','print','stylesheet','srinivas_latex.xsl');
%
%   publish2('Analysis_January.m',options);
% 

% close all to remove all extraneous figures, since we are publishing to a document anyway
close all