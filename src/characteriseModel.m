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
subplot(3,1,1), hold on
plot(time,data(1).response,'k')
l = plot(time,R,'r');
L = ['r^2=' oval(rsquare(data(1).response(a:z),R(a:z)))];
legend(l,L)
ylabel('Response (Hz)')
set(gca,'XLim',[a/1e3 z/1e3],'YLim',[min(data(1).response)-5 max(data(1).response)+5])
title(char(modelname))

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
end

subplot(3,3,4), hold on
ss  = 10;
for i = 1:length(data)
	plot(mean(data(i).stimulus(a:z)) + data(i).linear_prediction(a:ss:z),data(i).prediction(a:ss:z),'.','Color',c(i,:))
end
xlabel('Projected Stimulus (V)')
ylabel('Model Response (Hz)')

subplot(3,3,5), hold on
mean_stim = zeros(length(data),1);
for i = 1:length(data)
	plot(mean(data(i).stimulus(a:z)),data(i).gain,'+','Color',c(i,:))
	mean_stim(i) = mean(data(i).stimulus(a:z));
end

options = fitoptions(fittype('power1'));
options.Lower = [-Inf -1];
options.Upper = [Inf -1];
ff = fit(mean_stim,[data.gain]','power1',options);
plot(mean_stim,ff(mean_stim),'k--');

set(gca,'XScale','log','YScale','log')
ylabel('Gain (Hz/V)')
xlabel('Mean Stimulus (V)')
set(gca,'XLim',[min(min([data.stimulus]))*1.9 1.1*max(max([data.stimulus]))])

subplot(3,3,6), hold on
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

subplot(3,2,5), hold on
for i = length(contrasts):-1:1
	plot(contrast_data(i).linear_prediction(a:ss:z),contrast_data(i).resp(a:ss:z),'.','Color',c(i,:))
end
xlabel('Projected Stimulus (V)')
ylabel('Model Response (Hz)')

subplot(3,2,6), hold on
for i = length(contrasts):-1:1
	plot(std(contrast_data(i).stim),contrast_data(i).gain,'+','Color',c(i,:))
end
xlabel('Contrast')
ylabel('Gain (Hz/V)')
set(gca,'YLim',[.8*min([contrast_data.gain]) 1.2*max([contrast_data.gain])])


