% LVF_LFP.m
% recreates fig3 with LFP recordings
% 
% created by Srinivas Gorur-Shandilya at 11:52 , 14 July 2015. Contact me at http://srinivas.gs/contact/
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

%% 
% In this document, we attempt to recreate the results from Fig 3 showing fast gain changes with a widely distributed stimulus. 

%% Stimulus 
% First, we check that the stimulus is the same as the one we used six months ago (in LargeVarianceFlickering.m). The following plot shows the old stimulus (blue) and the new stimulus (red). 

% load old stimulus
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;

load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
PID = PID';
PID = PID(1:10:end,:);
% remove baseline
baseline = mean(PID(1:300,1));
PID = PID - baseline;

% load new stimulus
load('/local-data/DA-paper/large-variance-flicker/2015_07_14_LVFtest.mat')
newPID = data(4).PID';
newPID = newPID(1:10:end,:);
time = 1e-3*(1:length(PID));

% remove baseline
baseline = mean(newPID(1:300,1));
newPID = newPID - baseline;

% plot
clear l
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
l(1) = errorShade(time(1:10:end),mean2(PID(1:10:end,:)),sem(PID(1:10:end,:)),'Color',[0 0 1]);
l(2) = errorShade(time(1:10:end),mean2(newPID(1:10:end,:)),sem(newPID(1:10:end,:)),'Color',[1 0 0]);
xlabel('Time (s)')
ylabel('PID (V)')
set(gca,'XLim',[20 60])
legend(l,{'January 28 2015','July 14 2015'})

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% It looks very similar, but smaller, which is expected, as the PID gets less sensitive over time. What if we rescale it? 

ff = fit(mean2(newPID),mean2(PID),'poly1');

clear l
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
l(1) = errorShade(time(1:10:end),mean2(PID(1:10:end,:)),sem(PID(1:10:end,:)),'Color',[0 0 1]);
l(2) = errorShade(time(1:10:end),ff(mean2(newPID(1:10:end,:))),sem(newPID(1:10:end,:)),'Color',[1 0 0]);
xlabel('Time (s)')
ylabel('PID (V)')
legend(l,{'January 28 2015','July 14 2015 Rescaled'})
set(gca,'XLim',[20 60])
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

% tag the file as being published 
% add homebrew path
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1);

if being_published
	unix(strjoin({'tag -a published',which(mfilename)}));
end
