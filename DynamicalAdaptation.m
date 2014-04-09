%% Dynamical Adaptation in ORNs 
% Do ORNs exhibit fast adaptation to a flickering stimulus? Can a simple dynamical adaptation model predict ORN responses to flickering inputs? Here, I take data from Carlotta's random stimulus experiments and first check how well a simple linear prediction from the stimulus compares to the data. Then, I study how the instantaneous gain of the actual ORN response w.r.t to the predicted ORN response varies with time, and try to find a objective filter from the stimulus to this instantaneous gain to correct for systematic errors in the linear prediction.  

%% Adaptation vs. Response
% Most studies into the response and adaptation properties of sensory systems are conducted as follows: a long "conditioning" stimulus was presented to the sensor, and a "probe" stimulus was used to study its response. In Carlotta's experiments, she used a constant "background" odour concentration adn then applied pulses of odor on top of it and measured response to this. 
%
% <html>
% <img src="adaptation.PNG" width="800">
% </html>

%%
% Most neuron responses are generally non-linear to the amplitude of the stimulus applied. The LN model captures "gain compression" in the non-linear function of the model, and "true adaptaiton" is the differential response to pulses after these long conditioning pulses.

%%
% However, there is no fundamental basis for disambiguiating these two phenomenon. 

%% Adaptation and Response: Dynamical Adaptation
% Here, the adaptation and response of the neuron are considered together, and the neuron can be modeled as two filters that act on the input stimulus. One gives a linear response, and the other, the gain filter, modulates the linear response as a function of time in response to the stimulus. We expect the gain filter to be slower than the other, to account for adaptation. 
%
% <html>
% <img src="DA1.png" height="500">
% </html>

%% Synthetic Data: Sketch of Phenomenon
% In the following section, I will consider some synthetic data, and some synthetic filters, and see how these filters transform the inputs. I will then look to see if there is some characteristic relationships between the system output and the linear output. 

% parameters to tune for figure display
font_size = 20;
font_size2 = 18;
marker_size = 10;
marker_size2 = 20;

% define source
filename = '/data/random-stim/final_2011_06_14_ab3A_1o3ol3X-3_20ml_30sec_30ms_rand.mat';

% show input

%%
% Consider some hypothetical stimulus being presented to a ORN, which could look like this:
crop_traces = 0;
% grab data
if ~(exist('PID') == 1)
	crop_traces = 1;
	[PID, time, f, Valve] = PrepData3(filename);
	[PID1] = PrepData2(filename);

	% make all vectors consistent
	PID = PID(:);
	time = time(:);
	f = f(:);
	Valve = Valve(:);


	% correct offset
	PID = PID-min(PID);
	PID = PID/max(PID);
	PID = PID+min(PID1);
	PID = PID/max(PID);
	PID = PID*max(PID1);
	clear PID1

	% detrend PID
	ptrend = fit(time,PID,'Poly1'); 
	PID = PID - (ptrend(time) - mean(ptrend(time)));

	% detrend f
	ptrend = fit(time,f,'Poly1'); 
	f = f - (ptrend(time) - mean(ptrend(time)));

	clear ptrend

end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(time,PID,'r','LineWidth',3)
set(gca,'box','on','LineWidth',2,'XLim',[18 20],'FontSize',font_size)
ylabel('Stimulus Amplitude (a.u.)','FontSize',font_size)
xlabel('Time (s)','FontSize',font_size)
title('Synthetic Stimulus','FontSize',font_size)


% show filter
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


load smooth_fake_K.mat
% correct to ensure unit gain
fp = filter(K,1,PID-mean(PID)) + mean(f);
Kg = 0*K;
wl = 201;
Kg(33:33+wl) = 0.01;
gp = filter(Kg,1,PID-mean(PID)) + 1;
fpg = fp.*gp;
[ff] = fit(fpg,f,'Poly1');
K = K*ff.p1;
fp = filter(K,1,PID-mean(PID)) + mean(f);
fpg = fp.*gp;


%%
% Now, consider that the linear response of a ORN can be captured by a linear filter, that looks something like this:

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(filtertime,K,'LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[-0.11 0.91])
xlabel('Lag (s)','FontSize',font_size)
ylabel('Filter Amplitude','FontSize',font_size)
title('An example linear filter','FontSize',font_size)

% pass it through filter and show

%%
% Thus, the linear response of the ORN to this stimulus will look like this:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(time,fp,'k','LineWidth',3)
set(gca,'box','on','LineWidth',2,'XLim',[18 20])
ylabel('Filter Output (a.u.)','FontSize',font_size)


%%
% Now, consider that a second filter exists that modulates the gain of the previously described linear filter, and that it looks something like this: 

plot(filtertime,Kg,'LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[-0.11 0.91],'YLim',[-0.001 0.011])
xlabel('Lag (s)','FontSize',font_size)
ylabel('Filter Amplitude','FontSize',font_size)
title('Hypothetical gain filter','FontSize',font_size)

%%
% The output of the gain filter is the instantaneous modulation of the gain of the linear filter, and looks something like this:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(time,gp,'k','LineWidth',3)
set(gca,'box','on','LineWidth',2,'XLim',[18 20])
ylabel('Filter Output (a.u.)','FontSize',font_size)



%%
% The output of the linear filter (black) is compared to the gain modulated output (red). In some times, the output is lower, and in others, it is higher, depending on the average stimulus in the preceding window. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(time,fp,'k','LineWidth',3)
plot(time,fpg,'r','LineWidth',3)
set(gca,'box','on','LineWidth',2,'XLim',[18 20])
ylabel('Filter Output (a.u.)','FontSize',font_size)
legend LinearResponse GainModulated
title('Comparison of linear response and gain modulated output')

%%
% Now, the gain filter does not affect the linear filter at all, because the gain filter modulates a multiplicative factor on top of the linear filter output. We can check this by calculating a new filter from the PID to the output of the ORN, and comparing it to the linear filter, as follows:
%
% <html>
% <img src="DA2.png" height="500">
% </html>

K2 = FindBestFilter(PID,fpg);

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(filtertime,K,'g','LineWidth',2)
plot(filtertime,K2,'r','LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[-0.11 0.91])
xlabel('Lag (s)','FontSize',font_size)
ylabel('Filter Amplitude','FontSize',font_size)
title('Gain modulation does not affect filter','FontSize',font_size)

%% Synthetic Data: Signature of Dynamical Adaptation 
% In our toy model, where there are two filters, one controlling the gain of the other, we can characterise how the system responds to high and low stimuli as follows. We average the stimulus over some history window, and then look at how the system as a whole responds to the lowest 10% and highest 10% of smoothed stimuli, and compare it to the simple linear output of the system. 



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

history_lengths = [wl*3];
hl = history_lengths/3; % history lengths better be divisible by 3!
shat = NaN(length(hl),length(PID));
for i = 1:length(hl)
	shat(i,:) = filtfilt(ones(1,hl(i))/hl(i),1,PID);
	shat(i,1:hl(i)) = NaN;
end

[output_data] = GainAnalysis(fp,fpg,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,1);
xlabel('Linear Prediction','FontSize',font_size)
ylabel('Actual Model output (a.u.)','FontSize',font_size)


history_lengths = [30 102 150 201 300 402 501 600 801 1002];
hl = history_lengths/3; % history lengths better be divisible by 3!
shat = NaN(length(hl),length(PID));
for i = 1:length(hl)
	shat(i,:) = filtfilt(ones(1,hl(i))/hl(i),1,PID);
	shat(i,1:hl(i)) = NaN;
end

[output_data] = GainAnalysis(fp,fpg,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,[2 3]);


%% Real Data: Signature of Dynamical Adaptation 
% To check if this model is valid, or useful, we look for this signature of dynamical adaptation (higher gain to small stimuli than to large stimuli) in real data. To do this, we have to solve the following inverse problem:
% 
% <html>
% <img src="DA3.png" height="500">
% </html>

%% 
% Step 0: this is what we have to work with: the known stimulus to the ORN, and the known response from the ORN. 
figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
plot(time,PID,'k','LineWidth',3)
set(gca,'box','on','LineWidth',3,'XLim',[18 20],'FontSize',font_size)
ylabel('PID (a.u)','FontSize',font_size)
title('Real data: ab3A response to 1-octen-3-ol','FontSize',font_size)

subplot(2,1,2), hold on
plot(time,f,'k','LineWidth',3)
set(gca,'box','on','LineWidth',3,'XLim',[18 20],'FontSize',font_size)
ylabel('Firing Rate (Hz)','FontSize',font_size)



%%
% First, we find the linear filter from the stimulus to the ORN output. 
Kreal = FindBestFilter(PID,f);
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(filtertime,Kreal,'r','LineWidth',3)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[-0.11 0.91])
xlabel('Lag (s)','FontSize',font_size)
ylabel('Filter Amplitude','FontSize',font_size)
title('Filter: PID>ORN','FontSize',font_size)

%%
% Now that we know the linear filter, we can estimate the linear output using this filter. 
fp = filter(Kreal,1,PID-mean(PID)) + mean(f);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(time,f,'k','LineWidth',3)
plot(time,fp,'g','LineWidth',3)
set(gca,'box','on','LineWidth',2,'XLim',[18 20],'FontSize',font_size)
ylabel('Filter Output (a.u.)','FontSize',font_size)
xlabel('Time (s)','FontSize',font_size)
legend Data LinearPrediction

%%
% Now, we check if the same signature of dynamical adaptation is present in real data: 

history_lengths = [402];
hl = history_lengths/3; % history lengths better be divisible by 3!
shat = NaN(length(hl),length(PID));
for i = 1:length(hl)
	shat(i,:) = filtfilt(ones(1,hl(i))/hl(i),1,PID);
	shat(i,1:hl(i)) = NaN;
end


[output_data] = GainAnalysis(f,fp,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,1);

%%
% The following plot shows how the slope of the lines of best fit, or the instantaneous gains, varies with the history length. The plot on the right shows the goodness of fit for each fit, indicating the regions where the fit is meaningful.

history_lengths = [30 102 150 201 300 402 501 600 801 1002];
hl = history_lengths/3; % history lengths better be divisible by 3!
shat = NaN(length(hl),length(PID));
for i = 1:length(hl)
	shat(i,:) = filtfilt(ones(1,hl(i))/hl(i),1,PID);
	shat(i,1:hl(i)) = NaN;
end


[output_data] = GainAnalysis(f,fp,PID,shat,history_lengths,hl,filter_length,marker_size,marker_size2,font_size,[2 3]);

%% Real Data: Inference of Gain filter
% Now, we are given the input (the odour stimulus) and the output (the ORN response). We just calculated the linear filter and used it to estimate the linear output. The only things left are to find the gain filter and the gain modulation. 
% 
% <html>
% <img src="DA4.png" height="500">
% </html>

%%
% We can calculate the gain modulation from the output of the ORN and the linear prediction. However, that is noisy, and we can smooth it to remove high-frequency noise.  
gain2 = f./fp;
ws = [1 5 10 25 50 100];
gain2 = repmat(gain2,1,length(ws));
for i = 1:length(ws)
	gain2(:,i)=smooth(gain2(:,i),ws(i));	
end

figure('outerposition',[10 10 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold all
legendtext = {};
for i = 1:2:size(gain2,2)
	plot(time,gain2(:,i),'LineWidth',2)
	legendtext{i} = strcat('Window size:',mat2str(ws(i)*3),'ms');
end
set(gca,'XLim',[mean(time)-2 mean(time)+2],'LineWidth',2,'box','on','FontSize',font_size)
legendtext=legendtext(~cellfun('isempty',legendtext));
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

figure('outerposition',[10 10 500 500],'PaperUnits','points','PaperSize',[1200 500]); hold on
plot(ws*3,rvalues,'.','MarkerSize',marker_size2)
xlabel('Smoothing window (ms)','FontSize',font_size)
ylabel('r-square of corrected fit','FontSize',font_size)
line([1 max(ws)*3],[rsquare(f(filter_length+2:end),fp(filter_length+2:end)) rsquare(f(filter_length+2:end),fp(filter_length+2:end))],'LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size)

%%
% For with various smoothing of the instantaneous gain, we find the best gain filter for each and plot it below. Unfortunately, it looks pretty bad: 


if ~(exist('Kg_smooth') == 1)
	Kg_smooth = zeros(filter_length+1,size(gain2,2));
	for i = 1:size(gain2,2)
		Kg_smooth(:,i) = FindBestFilter(PID(filter_length+2:end),gain2(filter_length+2:end,i));
	end

end

figure('outerposition',[10 10 500 500],'PaperUnits','points','PaperSize',[1200 500]); hold on
plot(filtertime,Kg_smooth(:,1:2:end),'LineWidth',2), hold on
xlabel('Filter Lag','FontSize',font_size)
ylabel('Filter Amplitude','FontSize',font_size)
set(gca,'LineWidth',2,'FontSize',font_size,'box','on','XLim',[min(filtertime) max(filtertime)])
legend(legendtext)

%%
% How well can we predict the smoothed gain vectors using these filters and the stimulus? The blue dots show the theoretical maximum improvement of linear fit by gain modulation as a function of smoothing window. The black dots show the actual improvement (or lack thereof) from the estimation of the gain filter and the gain modulation. No point is above the horizontal line, which is the simple linear filter. 

gain2p = gain2;
for i = 1:length(ws)
	gain2p(:,i) = filter(Kg_smooth(:,i),1,PID-mean(PID)) + mean2(gain2(filter_length+2:end,i));
end

rvalues = NaN(1,length(ws));
for i = 1:length(ws)
	rvalues(i) = rsquare(gain2p(filter_length+2:end,i),gain2(filter_length+2:end,i));
end


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

figure('outerposition',[10 10 500 500],'PaperUnits','points','PaperSize',[1000 400]); hold on
plot(ws*3,rvalues,'b.','MarkerSize',marker_size2), hold on
plot(ws*3,rvalues2,'k.','MarkerSize',marker_size2)
xlabel('Smoothing window (ms)','FontSize',font_size)
ylabel('r-square of corrected fit','FontSize',font_size)
line([1 max(ws)*3],[rsquare(f(filter_length+2:end),fp(filter_length+2:end)) rsquare(f(filter_length+2:end),fp(filter_length+2:end))],'LineWidth',2)
set(gca,'LineWidth',2,'FontSize',font_size)


%% Real Data: Fitting the DA model 
% 
% <html>
% <img src="DAZ.png" height="500">
% </html>


close all