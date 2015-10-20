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
assert(strcmp(class(modelname),'function_handle'),'First argument must be a function handle to the model you want to characterise')

% first compare it to the lowest dose


R = modelname(data(1).stimulus,p);
time = 1e-3*(1:length(data(1).stimulus));

figure('outerposition',[0 0 900 1100],'PaperUnits','points','PaperSize',[900 1100]); hold on
subplot(4,1,1), hold on
plot(time,data(1).response,'k')
l = plot(time,R,'r');
L = ['r^2=' oval(rsquare(data(1).response(2e3:end),R(2e3:end)))];
legend(l,L)
ylabel('Response (Hz)')
set(gca,'XLim',[5 20])
title(char(modelname))

subplot(4,1,2:3), hold on
c = parula(length(data)+1);
for i = 1:length(data)
	data(i).prediction = modelname(data(i).stimulus,p);
	plot(time,data(i).prediction,'Color',c(i,:))
end
set(gca,'XLim',[5 20])
xlabel('Time (s)')
ylabel('Response (Hz)')

% extract filters in all cases
for i = 1:length(data)
	[K(:,i), prediction, gain, gain_err] = extractFilters(data(i).stimulus,data(i).prediction,'a',2e3,'z',18e3);
	data(i).linear_prediction = prediction;
	data(i).gain = gain;
	gata(i).gain_err = gain_err;
end

subplot(4,3,10), hold on
ss  = 10;
for i = 1:length(data)
	plot(mean(data(i).stimulus(2e3:end)) + data(i).linear_prediction(2e3:ss:end),data(i).prediction(2e3:ss:end),'.','Color',c(i,:))
end
xlabel('Projected Stimulus (V)')
ylabel('Model Response (Hz)')

subplot(4,3,11), hold on
for i = 1:length(data)
	plot(mean(data(i).stimulus(2e3:end)),data(i).gain,'+','Color',c(i,:))
end

ff = fit(mean([data.stimulus])',[data.gain]','power1');
l = plot(mean([data.stimulus]),ff(mean([data.stimulus])),'k--');
L = ['\alpha = ' oval(ff.b)];
legend(l,L)
set(gca,'XScale','log','YScale','log')
ylabel('Gain (Hz/V)')
xlabel('Mean Stimulus (V)')

subplot(4,3,12), hold on
filtertime = (1:length(K)) - 100;
for i = 1:length(data)
	plot(filtertime,K(:,i),'Color',c(i,:))
end
xlabel('Filter Lag (ms)')



