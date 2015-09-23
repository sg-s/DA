% SwitchingAnalysis.m
% This document analyses stimuli to be used in the switching experiment, where we switch between one variance/mean and another over and over again.
% 
% created by Srinivas Gorur-Shandilya at 1:20 , 25 August 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

% this code determines if this function is being called by publish() or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic

%% Switching Analysis.m
% This document analyses stimuli to be used in the switching experiment, where we switch between one variance/mean and another over and over again.

%% MFC configuration
% In this document, we configure the MFC so that it runs the PD algorithm, with parameters P = 2500 and D = 10000. With these values, the MFC is very fast (approx. 20ms response timescale), but not reproducible. There is noticeable ringing on step changes, but they die out quickly. 

%%
% Control signals are generated using MakeMeanSwitching.m and MakeVarianceSwitching.m. 

load('/local-data/DA-paper/switching/2015_08_25_MeanSwitchingTest.mat')

%% Mean Switching 
% We now look at the stimuli for the mean switching case. In the following figure, we look at the raw trace, to get a sense of what the data looks like. 

figure('outerposition',[0 0 900 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
subplot(3,1,1),hold on
t = 1e5:10:2.5e5;
plot(1e-4*t,ControlParadigm(1).Outputs(1,t),'k')
ylabel('MFC Setpoint (V)')
subplot(3,1,2),hold on
plot(1e-4*t,data(1).MFC500(1,t),'k')
ylabel('MFC Flow (V)')
subplot(3,1,3),hold on
plot(1e-4*t,data(1).PID(1,t),'k')
ylabel('PID (V)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%% 
% As we can see from this example trace, even though the MFC is doing a good job keeping up with the command signal, the PID values change slowly. We now make the same plot, but averaged over all the switches. 

a = 5e4+1:10e4:length(data(1).PID);
cs = NaN(15e3,length(a));
mfc = cs;
pid = cs;

for i = 2:length(a)
	this_start = a(i);
	try
		temp = ControlParadigm(1).Outputs(1,this_start-5e4:this_start+10e4);
		cs(:,i) = temp(1:10:end-1);

		temp = data(1).MFC500(1,this_start-5e4:this_start+10e4);
		mfc(:,i) = temp(1:10:end-1);

		temp = data(1).PID(1,this_start-5e4:this_start+10e4);
		pid(:,i) = temp(1:10:end-1);
	end
end

figure('outerposition',[0 0 900 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
subplot(3,1,1),hold on
errorShade(1e-3*(1:length(cs)),mean2(cs),sem(cs),'Color',[0 0 0],'Shading',0.7);
ylabel('MFC Setpoint (V)')
subplot(3,1,2),hold on
errorShade(1e-3*(1:length(mfc)),mean2(mfc),sem(mfc),'Color',[0 0 0],'Shading',0.7);
ylabel('MFC Flow (V)')
subplot(3,1,3),hold on
errorShade(1e-3*(1:length(pid)),mean2(pid),sem(pid),'Color',[0 0 0],'Shading',0.7);
ylabel('PID (V)')

PrettyFig()

load('/local-data/DA-paper/switching/2015_08_25_VarianceSwitchingTest.mat')

%% Variance Switching 
% We now look at the stimuli for the variance switching case. In the following figure, we look at the raw trace, to get a sense of what the data looks like. 

figure('outerposition',[0 0 900 900],'PaperUnits','points','PaperSize',[1200 900]); hold on
subplot(3,1,1),hold on
t = 1e5:10:2.5e5;
plot(1e-4*t,ControlParadigm(1).Outputs(1,t),'k')
ylabel('MFC Setpoint (V)')
subplot(3,1,2),hold on
plot(1e-4*t,data(1).MFC500(1,t),'k')
ylabel('MFC Flow (V)')
subplot(3,1,3),hold on
plot(1e-4*t,data(1).PID(1,t),'k')
ylabel('PID (V)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% We now plot the distributions of the signals in the two regimes:

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
x = 0:0.01:1.5;
y = NaN(length(x),20);
a = 5e4+1:10e4:length(data(1).PID);
for i = 1:length(a)
	this_start = a(i);
	try
		temp = ControlParadigm(1).Outputs(1,this_start-5e4:this_start);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[0 0 1],'Shading',0.4);
y = NaN(length(x),20);
a = 10e4+1:10e4:length(data(1).PID);
for i = 1:length(a)
	this_start = a(i);
	try
		temp = ControlParadigm(1).Outputs(1,this_start-5e4:this_start);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[1 0 0],'Shading',0.4);
xlabel('MFC Setpoint (V)')
ylabel('Probability')

subplot(1,3,2), hold on
x = 0:0.01:1.5;
y = NaN(length(x),20);
a = 5e4+1:10e4:length(data(1).PID);
for i = 1:length(a)
	this_start = a(i);
	try
		temp = data(1).MFC500(1,this_start-5e4:this_start);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[0 0 1],'Shading',0.4);
y = NaN(length(x),20);
a = 10e4+1:10e4:length(data(1).PID);
for i = 1:length(a)
	this_start = a(i);
	try
		temp = data(1).MFC500(1,this_start-5e4:this_start);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[1 0 0],'Shading',0.4);
xlabel('MFC Flow (V)')

subplot(1,3,3), hold on
x = 0:0.01:2;
y = NaN(length(x),20);
a = 5e4+1:10e4:length(data(1).PID);
for i = 1:length(a)
	this_start = a(i);
	try
		temp = data(1).PID(1,this_start-5e4:this_start);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[0 0 1],'Shading',0.4);
y = NaN(length(x),20);
a = 10e4+1:10e4:length(data(1).PID);
for i = 1:length(a)
	this_start = a(i);
	try
		temp = data(1).PID(1,this_start-5e4:this_start);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[1 0 0],'Shading',0.4);
xlabel('Stimulus (V)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end

t = toc;

%% 
% This document was built in: 
disp(strcat(oval(t,3),' seconds.'))

% tag the file as being published 
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(strjoin({'tag -a published',which(mfilename)}));
end
