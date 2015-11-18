% characteriseModel.m
% accepts a model, some best-fit parameters, and characterises how this model responds to stimuli of increasing mean
% meant to by called by FittingModels2Data.m
% 
% created by Srinivas Gorur-Shandilya at 10:57 , 20 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = characteriseModel(modelname,p,data)

% defensive programming 
assert(isa(modelname,'function_handle'),'First argument must be a function handle to the model you want to characterise')

% first compare it to the lowest dose
R = modelname(data(1).stimulus,p);
time = 1e-3*(1:length(data(1).stimulus));

a = 5e3;
z = length(data(1).response);

figure('outerposition',[0 0 900 1100],'PaperUnits','points','PaperSize',[900 1100]); hold on
suptitle(char(modelname))

% generate predictions for all means
c = parula(length(data)+1);
for i = 1:length(data)
	data(i).prediction = modelname(data(i).stimulus,p);
end

% extract filters in all cases
K = zeros(700,length(data));
for i = 1:length(data)
	[K(:,i), prediction, gain, gain_err] = extractFilters(data(i).stimulus,data(i).prediction,'a',a,'z',z);
	data(i).linear_prediction = prediction;
	data(i).gain = gain;
	data(i).gain_err = gain_err;
	data(i).gain2 = std(data(i).prediction(a:z))/std(data(i).stimulus(a:z));
end

subplot(3,3,1), hold on
ss  = 10;
for i = 1:length(data)
	plot(mean(data(i).stimulus(a:z)) + data(i).linear_prediction(a:ss:z),data(i).prediction(a:ss:z),'.','Color',c(i,:))
end
xlabel('Projected Stimulus (V)')
ylabel('Model Response (Hz)')

subplot(3,3,2), hold on
mean_stim = zeros(length(data),1);
for i = 1:length(data)
	plot(mean(data(i).stimulus(a:z)),data(i).gain,'+','Color',c(i,:))
	mean_stim(i) = mean(data(i).stimulus(a:z));
end

options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
gain = [data.gain];
mean_stim(isnan(gain)) = [];
gain(isnan(gain)) = [];
ff = fit(mean_stim(:),gain(:),'power1',options);
plot(mean_stim,ff(mean_stim),'k--');

set(gca,'XScale','log','YScale','log')
ylabel('Gain (Hz/V)')
xlabel('Mean Stimulus (V)')
set(gca,'XLim',[min(min([data.stimulus]))*1.9 1.1*max(max([data.stimulus]))])

subplot(3,3,3), hold on
[~,loc]=max(K);
plot(mean([data.stimulus]),loc-100,'k+')
ylabel('Response Timescale (ms)')
xlabel('Mean Stimulus (V)')
set(gca,'YLim',[0 1.1*max(loc-100)])

% generate responses do different contrasts
contrasts = [.1 .2 .3 .5 1 1.2 1.5 2];
clear contrast_data
contrast_data.gain = NaN;
contrast_data.resp = []; 
contrast_data.linear_prediction = [];
K2 = zeros(700,length(contrasts));
c = parula(length(contrasts)+1);
for i = 1:length(contrasts)
	contrast_data(i).stim = (1+(data(3).stimulus - mean(data(3).stimulus))*contrasts(i));
	contrast_data(i).resp = modelname(contrast_data(i).stim,p);
	[K2(:,i), linear_prediction, gain] = extractFilters(contrast_data(i).stim,contrast_data(i).resp,'a',a,'z',z);
	contrast_data(i).linear_prediction = linear_prediction;
	contrast_data(i).gain = gain;
end

subplot(3,2,3), hold on
for i = length(contrasts):-1:1
	plot(contrast_data(i).linear_prediction(a:ss:z),contrast_data(i).resp(a:ss:z),'.','Color',c(i,:))
end
xlabel('Projected Stimulus (V)')
ylabel('Model Response (Hz)')

subplot(3,2,4), hold on
for i = length(contrasts):-1:1
	plot(std(contrast_data(i).stim),contrast_data(i).gain,'+','Color',c(i,:))
end
xlabel('Contrast')
ylabel('Gain (Hz/V)')
set(gca,'YLim',[.8*min([contrast_data.gain]) 1.2*max([contrast_data.gain])])


% generate triangles 
Tmax = 20;
peak_locs = 1:2:Tmax;
c = parula(length(peak_locs) + 1);
ax(1) = subplot(3,2,5); hold(ax(1),'on')
title('Triangle Stimulus')
xlabel('Time (s)')
ax(2) = subplot(3,2,6); hold(ax(2),'on')
title('Model Response')
ylabel('Response (Hz)')
xlabel('Time (s)')
t = 1e-3:1e-3:Tmax;
for i = 1:length(peak_locs)
	s = 0*t;
	rs = 1e-3/peak_locs(i);
	s(1:peak_locs(i)*1e3) = rs:rs:1;
	rs = 1e-3/(Tmax-peak_locs(i));
	s(peak_locs(i)*1e3+1:end) = 1:-rs:rs;
	plot(ax(1),t,s,'Color',c(i,:))

	r = modelname(s,p);
	plot(ax(2),t,r,'Color',c(i,:))
end


