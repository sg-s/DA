% Flickering Stimulus on a Background
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

if ~exist('data')
	load('/local-data/DA-paper/flickering-stim-back/pid/2014_09_23_srinivas_2ac_flickering_background.mat')
end

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end

% remember to set redo_bootstrap

% ########  #### ########  
% ##     ##  ##  ##     ## 
% ##     ##  ##  ##     ## 
% ########   ##  ##     ## 
% ##         ##  ##     ## 
% ##         ##  ##     ## 
% ##        #### ########  

%%
% Carlotta's data is almost entirely in the liner regime, and the effects of fast gain control, if any, are very small. Mahmut's "natural" stimuli are too poorly controlled, and the neuron response can't be fit by any model we have, making analysis challenging.

%%
% This document shows a middle path: we use a random flickering stimulus as in Carlotta's experiments, but now the flickering stimuli have varying heights, in addition to varying durations. Furthermore, they occur on top of a background, whose height can be varied too. 

%%
% The following figure shows one such paradigm, where pulses are presented on top of a background. The odor used is Ethyl Acetate. The different colors are different trials of the same experiment, and the black line is a fit to the background, showing good trial-to-trial reproducibility. 


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
time = 1e-4:1e-4:1e-4*length(data(1).PID);
c = jet(width(data(1).PID));
for i = 1:width(data(1).PID)
	plot(time,data(1).PID(i,:),'Color',c(i,:))
end
set(gca,'XLim',[4 14])

% figure out the background level and show it
valve = ControlParadigm(1).Outputs(5,:);
[ons,offs]=ComputeOnsOffs(valve);
% for each segment, find the minimum
minPID = NaN*ons;
t_min = NaN*ons;
PID = mean2(data(1).PID);
for i = 1:length(ons)
	[minPID(i),loc]=min(PID(ons(i):offs(i)));
	t_min(i) = loc+ons(i);
end

c =mean(minPID)+.5*std(minPID);
t_min(minPID>c)=[];
minPID(minPID>c)=[];
% fit a line
ff = fit(t_min*1e-4,minPID,'poly1');
plot(time,ff(time),'k')
xlabel('Time (s)')
ylabel('PID (V)')

title(strcat('m=',oval(ff.p1,2),'V/s'))

PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%% 
% The following figure shows the same stimulus on top of backgrounds of different amplitudes. Once again, we see that the background is mostly constant, and the pulse height is not really affected by the background. The lines are best fits to the estimated background concentrations. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = jet(length(data));
for i = 1:length(data)
	PID = mean(data(i).PID);
	plot(time,PID,'Color',c(i,:))

	% figure out the background level and show it
	valve = ControlParadigm(i).Outputs(5,:);
	[ons,offs]=ComputeOnsOffs(valve);
	% for each segment, find the minimum
	minPID = NaN*ons;
	t_min = NaN*ons;
	for j = 1:length(ons)
		[minPID(j),loc]=min(PID(ons(j):offs(j)));
		t_min(j) = loc+ons(j);
	end

	cutoff =mean(minPID)+.5*std(minPID);
	t_min(minPID>cutoff)=[];
	minPID(minPID>cutoff)=[];
	% fit a line
	ff = fit(t_min*1e-4,minPID,'poly1');
	plot(time,ff(time),'Color',c(i,:))

end

xlabel('Time (s)')
ylabel('PID (V)')
title('Ethyl acetate pulses on top of backgrounds')
set(gca,'XLim',[12 18],'box','on')
PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end

%% Synthetic Data
% In this section, we generate a synthetic neuron response using a DA model with parameters chosen to exhibit a large "gain adaptation". We then fit a LN model to this and determine if this experimental paradigm is right for observing the effects of gain control, if any, in this dataset. The parameters of the DA model were chosen to match the dose-response properties of ab3A. The parameters are:

% parameters of DA model
p.s0 = 0;
p.tau_y = 1; % assuming dt=3e-3
p.tau_z = 7;
p.n_y = 2;
p.n_z = 2;
p.C = 0.3;
p.A = 300;
p.B = 3;

disp(p)

for i = 1:length(data)
	t = 3e-3:3e-3:max(time);
	thisPID = interp1(time,mean2(data(i).PID),t);
	p.s0  = min(thisPID(600:1600));
	data(i).DAModel = DA_integrate2(thisPID,p);
end

%%
% The following figure shows the response of the DA model with the chosen parameters to the stimulus sets:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = jet(length(data));
for i = 1:length(data)
	plot(t,data(i).DAModel,'Color',c(i,:))
end
xlabel('Time (s)')
ylabel('Response (Hz)')
title('Response of DA Model to stimulus')
set(gca,'XLim',[12 18],'YLim',[0 200],'box','on')
PrettyFig;

if being_published
	snapnow;
	delete(gcf);
end


%%
% We now fit a LN model to each of these traces. The following plot shows the best-fit LN model predictions to two DA Model simulations (for the lowest background and the highest background). 


a = 3000;
z = 8000;
for i = 1:length(data)
	t = 3e-3:3e-3:max(time);
	thisPID = interp1(time,mean2(data(i).PID),t);

	% build a simple linear model
	[data(i).K,~,filtertime] = FindBestFilter(thisPID(a:z),data(i).DAModel(a:z),[],'filter_length=201;');
	data(i).filtertime = filtertime*mean(diff(t));
	data(i).LinearFit = convolve(t,thisPID,data(i).K,data(i).filtertime);
	data(i).LinearFit = data(i).LinearFit + mean(data(i).DAModel(a:z));

	xdata = data(i).LinearFit(a:z);
	ydata = data(i).DAModel(a:z);

	% crop it to lose NaNs
	ydata(isnan(xdata)) = [];
	xdata(isnan(xdata)) = [];

	xdata = xdata(:);
	ydata = ydata(:);

	fo=optimset('MaxFunEvals',1000,'Display','none');
	x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);
	% save this for later
	data(i).LN_pred = hill(x,data(i).LinearFit);
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = [1 length(data)]
	plot(t,data(i).DAModel,'Color','k')
	plot(t,data(i).LN_pred,'Color','r')
end
legend({'Simulated Data','LN Fit'})
xlabel('Time (s)')
ylabel('Response (Hz)')
title('Performance of LN Model to simulated responses on background')
set(gca,'XLim',[12 18],'YLim',[0 200],'box','on')
PrettyFig;
if being_published
	snapnow;
	delete(gcf);
end

%%
% Under what regimes does the LN model do well? In the following figure, we plot the success of the LN model (as r-square) vs. the background odor level. 

background_odor_level = NaN(length(data),1);
r = NaN(length(data),1);
for i = 1:length(data)
	thisPID = interp1(time,mean2(data(i).PID),t);
	background_odor_level(i) = min(thisPID(a:z));
	r(i) = rsquare(data(i).DAModel(a:z),data(i).LN_pred(a:z));
end
clear i

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
scatter(background_odor_level,r)
xlabel('Background Level (V)')
ylabel('LN Model Fit Quality (r-square)')
PrettyFig;
if being_published
	snapnow;
	delete(gcf);
end

%% In which regimes do LN model predictions hint gain adaptation? 
% Poor LN Model fits are indicative of interesting dynamics that cannot be explained by models that excel in stationary stimuli (like a LN Model). On the other hand, we know that LN models are very good when there is no background, and when the neuron silences between pulses. Thus, the stimulus we want to deliver is the one with the smallest background such that the neuron does not go silent between pulses. 

%%
% In the following figures, we perform a gain analysis using standard methods described elsewhere in this project on the LN model predictions of the simulated responses to the stimulus with the smallest background. 

history_lengths = (3*floor(1000*logspace(-2,1,30)/3))/1e3;
example_history_length = 0.135;
td=1;
f1=figure('outerposition',[0 0 1000 600],'PaperUnits','points','PaperSize',[1000 600]); hold on
ph(1) = subplot(2,1,1); hold on 
ph(2) = subplot(2,1,2); hold on

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on


s = 3000; % when we start for the gain analysis
z = 8000; % where we end

clear x
x.response = data(td).DAModel(s:z);
x.prediction = data(td).LN_pred(s:z);
thisPID = interp1(time,mean2(data(td).PID),t);
x.stimulus = thisPID(s:z);
x.time = t(s:z);
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


%%
% How do the measured effects of gain adaptation vary with levels of background stimuli? The following figure shows various relative gain measurements for the different background levels: 

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
if redo_bootstrap
	all_p = zeros(length(history_lengths),length(data));
end
yy=[];
yy(1) = Inf; yy(2) = -Inf;
for td = 1:length(data)
	ph = [];
	ph(4) = subplot(2,3,td); hold on

	clear x
	x.response = data(td).DAModel(s:z);
	x.prediction = data(td).LN_pred(s:z);
	thisPID = interp1(time,mean2(data(td).PID),t);
	x.stimulus = thisPID(s:z);
	x.time = t(s:z);
	x.filter_length = 201;

	if redo_bootstrap
		temp = GainAnalysis4(x,history_lengths,example_history_length,ph);
		temp = temp(1,:);
		all_p(:,td) = temp;

	else
		GainAnalysis4(x,history_lengths,example_history_length_LN,ph, [all_p(:,td) all_p(:,td)]');
	end
	set(ph(4),'XScale','log')
	this_yy = get(ph(4),'YLim');
	yy(1) = min([yy(1) this_yy(1)]);
	yy(2) = max([yy(2) this_yy(2)]);
end

for td = 1:length(data)
	subplot(2,3,td)
	set(gca,'YLim',yy)
	titlestr = strcat('Background :',oval(background_odor_level(td),2));
	title(titlestr)
end

%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))



