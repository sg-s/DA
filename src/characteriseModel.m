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
subplot(4,1,1), hold on
plot(time,data(1).response,'k')
l = plot(time,R,'r');
L = ['r^2=' oval(rsquare(data(1).response(a:z),R(a:z)))];
legend(l,L)
ylabel('Response (Hz)')
set(gca,'XLim',[a/1e3 z/1e3],'YLim',[min(data(1).response)-5 max(data(1).response)+5])
title(char(modelname))

subplot(4,1,2:3), hold on
c = parula(length(data)+1);
for i = 1:length(data)
	data(i).prediction = modelname(data(i).stimulus,p);
	plot(time(a:z),data(i).prediction(a:z),'Color',c(i,:))
end
set(gca,'XLim',[a/1e3 z/1e3])
xlabel('Time (s)')
ylabel('Response (Hz)')

% extract filters in all cases
K = zeros(700,length(data));
for i = 1:length(data)
	[K(:,i), prediction, gain, gain_err] = extractFilters(data(i).stimulus,data(i).prediction,'a',a,'z',z);
	data(i).linear_prediction = prediction;
	data(i).gain = gain;
	data(i).gain_err = gain_err;
end

subplot(4,3,10), hold on
ss  = 10;
for i = 1:length(data)
	plot(mean(data(i).stimulus(a:z)) + data(i).linear_prediction(a:ss:z),data(i).prediction(a:ss:z),'.','Color',c(i,:))
end
xlabel('Projected Stimulus (V)')
ylabel('Model Response (Hz)')

subplot(4,3,11), hold on
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

subplot(4,3,12), hold on
[~,loc]=max(K);
plot(mean([data.stimulus]),loc-100,'k+')
ylabel('Response Timescale (ms)')
xlabel('Mean Stimulus (V)')




