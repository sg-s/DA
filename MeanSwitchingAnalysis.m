% MeanSwitchingAnalysis.m
% analysis of mean switching experiment
% 
% created by Srinivas Gorur-Shandilya at 2:34 , 27 August 2015. Contact me at http://srinivas.gs/contact/
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

%% Mean Switching Analysis
% In this document, we analyse the responses of ORNs to stimuli where we rapidly switch from two ensembles of stimuli, differing only in their means. 

%% Stimulus
% In the following figure, we show what the stimulus looks like: 

p = '/local-data/switching/mean/use-this/';
[PID, LFP, fA, paradigm, orn] = consolidateData(p,1);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,4,1:3), hold on
time = 1e-3*(1:length(PID));
plot(time,PID(:,1),'k');
set(gca,'XLim',[30 100])
xlabel('Time (s)')
ylabel('Stimulus (V)')

subplot(1,4,4), hold on
x = 0:0.01:2.5;
y = NaN(length(x),width(PID));
a = 10e3+1:10e3:length(PID)-5e3;
for i = 1:length(a)
	this_start = a(i);
	try
		temp = PID(this_start-5e3:this_start,2);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[1 0 0],'Shading',0.4);
x = 0:0.01:2.5;
y = NaN(length(x),width(PID));
a = 15e3+1:10e3:length(PID)-5e3;
for i = 1:length(a)
	this_start = a(i);
	try
		temp = PID(this_start-5e3:this_start,2);
		y(:,i) = hist(temp(1:10:end),x);
		y(:,i) = y(:,i)/sum(y(:,i));
	end
end
errorShade(x,mean2(y),sem(y),'Color',[0 0 1],'Shading',0.4);
xlabel('Stimulus (V)')
ylabel('Probability')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end


%% LFP Analysis
% We now line up all the LFP signals by the switching time, and look at how the LFP changes in aggregate when we switch from a low mean stimulus to a high mean stimulus. Since the LFP drifts over time, and we don't really care about that, we construct a triggered LFP where we subtract the mean of the LFP during the low stimulus window in each epoch and then average. 


triggeredLFP = NaN(1e4,0);
a = 20e3;
z = length(LFP)-5e3;
for start_here = a:10e3:z
	for i = 1:width(LFP)
		this_segment = LFP(start_here:start_here+1e4-1,i);
		this_segment = this_segment - mean(this_segment(1:5e3));
		triggeredLFP = [triggeredLFP this_segment];
	end
end

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorShade(1e-3*(1:length(triggeredLFP)),mean2(triggeredLFP),sem(triggeredLFP))
xlabel('Time since high \rightarrow low switch (s)')
ylabel('\DeltaLFP (mV)')

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
