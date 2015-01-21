% Paper Figures
% makes all the figures for the paper
% 
% created by Srinivas Gorur-Shandilya at 12:57 , 21 January 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
% redo = 1; deliberately unset

% intially we will use the trace from the mean shufted gaussian experiment
load('MeanShiftedGaussians.mat')
plot_these=find(strcmp(paradigm_names{1}, combined_data.paradigm));
% remove trend
b = floor(5/3e-3);
a = floor(15/3e-3);
z = floor(55/3e-3);
clear detrended_data
detrended_data.time = [];
detrended_data.stim = [];
detrended_data.resp = [];
for i = 1
	plot_these=find(strcmp(paradigm_names{i}, combined_data.paradigm));
	this_pid=mean2(combined_data.PID(plot_these,:));
	this_resp=mean2(combined_data.fA(:,plot_these));
	time = 3e-3*(1:length(this_resp));
	baseline = mean(this_pid(1:b));
	time = time(a:z);
	this_pid = this_pid(a:z);
	this_resp = this_resp(a:z);

	detrended_data(i).time = time;
	ff = fit(time(:),this_pid(:),'poly2');
	detrended_data(i).stim = this_pid - ff(time)' + mean(ff(time)) - baseline;

	ff = fit(time(:),this_resp(:),'poly2');
	detrended_data(i).resp = this_resp - ff(time) + mean(ff(time));
end
dt = mean(diff(detrended_data.time));

%% Figure 1
% Olfactory Receptor Neurons control gain on a fast time scale in response to fluctuating odor stimuli. 
figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
a(1) = subplot(9,3,[1 2 4 5]); hold on % stimulus
a(2) = subplot(9,3,[7 8 10 11]); hold on % response + prediction
a(3) = subplot(9,3,[13 14 16 17]); hold on% gain

% LN Model
a(4) = subplot(9,3,[3 6 9]);hold on
a(5) = subplot(9,3,[12 15 18]);hold on

% gain analysis
a(6) = subplot(9,3,[19:3:25]);hold on % scatter plot
a(7) = subplot(9,3,[20:3:26]); hold on% history lengths plot
a(7) = subplot(9,3,[21:3:27]); hold on% autocorrelation

% make the plots
plot(a(1),detrended_data.time,detrended_data.stim,'k')
ylabel(a(1),'Stimulus (V)')
plot(a(2),detrended_data.time,detrended_data.resp,'k')
ylabel(a(2),'ORN Response (Hz)')

p.  tau1= 6.0236;
p.   K_n= 3.6445;
p.  tau2= 22.7344;
p.   K_A= 0.4189;
p.     A= 45.6247;
p.     n= 3.1602;
p.    Kd= 16.5625;
p.offset= 15.8848;

fp = pLNModel(detrended_data.stim,p);
plot(a(2),detrended_data.time,fp,'r')

plot(a(3),detrended_data.time,detrended_data.resp./fp','k')
ylabel(a(3),'Gain')

% plot the LN model
t = 1:300;
K = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);
t = dt*t;
plot(a(4),t,K,'k')
xlabel(a(4),'Filter Lag (s)')

x = [p.A p.Kd p.n];
f = hill(x,[0:50]);
plot(a(5),[0:50],f,'k')




%% Figure 2
% Fast Gain Control is a general phenomenon, seen in different receptor-odor combinations 