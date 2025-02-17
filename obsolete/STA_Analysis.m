% STA_Analysis.m
% 
% created by Srinivas Gorur-Shandilya at 4:18 , 18 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end
tic


%% STA Analysis
% In this document we try to perform a spike-triggered analysis of the data to see if we can quantify the data without resorting to firing rates, which inherently smooth and average data. 

%%
% In the following figures, we plot the spike triggered average for two sensilla. The odour is ethyl acetate, which is presented in a flickering fashion (the data from the large variance flicker is used here). Each figure corresponds to one sensilla, and we show the analysis for both A and B neurons. In this analysis, we normalise the stimulus by subtracting the mean and dividing through by the standard deviation. 

before = 1e4;
after = 6e3;

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
subplot(2,1,1), hold on
K = STA(spikes(4).A,data(4).PID,before,after);
t = -before:after;
t = t*1e-4;
plot(t,K)
ylabel('Stimulus (norm)')
title('A neuron')

subplot(2,1,2), hold on

K = STA(spikes(4).B,data(4).PID,before,after);
t = -before:after;
t = t*1e-4;
plot(t,K)
ylabel('Stimulus (norm)')
title('B neuron')
xlabel('Time relative to spike (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
K = STA(spikes(4).A,data(4).PID,before,after);
t = -before:after;
t = t*1e-4;
plot(t,K)
ylabel('Stimulus (norm)')
title('A neuron')

subplot(2,1,2), hold on

K = STA(spikes(4).B,data(4).PID,before,after);
t = -before:after;
t = t*1e-4;
plot(t,K)
ylabel('Stimulus (norm)')
title('B neuron')
xlabel('Time relative to spike (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% We get very clean filters without having to calculate firing rates at all. Now, we repeat the analysis on Carlotta's flickering data. 

before = 1e4;
after = 6e3;
load('/local-data/DA-paper/carlotta-martelli/flickering-stim/data.mat')
figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
K = STA(data(1).spikes(1e5:end-5e4,:),data(1).fullPID(1e5:end-5e4,:),before,after);
t = -before:after;
t = t*1e-4;
plot(t,K)
title(strrep(data(1).original_name,'_','-'))
ylabel('Stimulus (norm)')

subplot(2,1,2), hold on
K = STA(data(8).spikes(1e5:end-2e5,:),data(8).fullPID(1e5:end-2e5,:),before,after);
t = -before:after;
t = t*1e-4;
plot(t,K)
title(strrep(data(8).original_name,'_','-'))
ylabel('Stimulus (norm)')
xlabel('Time relative to spike (s)')

PrettyFig;

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


