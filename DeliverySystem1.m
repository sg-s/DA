% DeliverySystem1.m
% First attempt at understanding how the delivery system works
% created by Srinivas Gorur-Shandilya at 10:21 , 22 July 2015. Contact me at http://srinivas.gs/contact/
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


%% Turbulence and PID position
% In this document we look at effects of turbulence and PID position and how that affects our ability to measure the stimulus. 

%%
% The experimental setup is shown in the figure below. Breifly, a glass tube with three holes is used. The most upstream hole has a valve, and the PID is inserted into the most downstream hole. All the holes are on the same side. In some experiments, a coiled wire is inserted at the mouth of the tube (upstream) in an attempt to introduce turbulence into the airstream so that it mixes more. 
% 
% <</Users/sigbhu/code/da/images/2015-7-22-fig1.png>>
%

%%
% In the following figure, we show the flow rate distributions and stimulus distributions from the three positions of the PID, for the case where there is no mixing wire: 

load('/local-data/DA-paper/MSG-construction/2015_07_22_exp1_sans_mixer.mat')
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
pos = [3 3 1 1 2 2];
c = parula(4);
for i = 1:width(data.MFC500)
	[y,x] = hist(data.MFC500(i,20e4:10:55e4),50);
	y = y/sum(y);
	plot(x*100,y,'Color',c(pos(i),:))
end
xlabel('Flow Rate (mL/min)')
subplot(1,2,2), hold on
c = parula(4);
clear l
for i = 1:width(data.PID)
	[y,x] = hist(data.PID(i,20e4:10:55e4),50);
	y = y/sum(y);
	l(pos(i)) = plot(x,y,'Color',c(pos(i),:));
end
xlabel('Stimulus (V)')
legend(l,{'Pos1','Pos2','Pos3'},'location','northeast')
PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

olddata = data;

%% 
% Now we introduce the mixing wire and see what effect this has:


load('/local-data/DA-paper/MSG-construction/2015_07_22_exp1.mat')
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
pos = [1 1 3 3 2 2];
c = parula(4);
clear l
for i = 1:width(data.MFC500)
	[y,x] = hist(data.MFC500(i,20e4:10:55e4),50);
	y = y/sum(y);
	plot(x*100,y,'Color',c(pos(i),:))
end
xlabel('Flow Rate (mL/min)')
subplot(1,2,2), hold on
c = parula(4);
for i = 1:width(data.PID)
	[y,x] = hist(data.PID(i,20e4:10:55e4),50);
	y = y/sum(y);
	l(pos(i)) = plot(x,y,'Color',c(pos(i),:));
end
xlabel('Stimulus (V)')
legend(l,{'Pos1','Pos2','Pos3'},'location','northeast');
PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% These experiments stress the importance of mixing the secondary air stream nicely into the main air stream, and also of having a precisely defined PID position, since this clearly makes a huge difference. 

%% Mixer below valve
% We now attempt to mix the streams better by having the mixer below the valve, and upstream of the PID

load('/local-data/DA-paper/MSG-construction/2015_07_22_exp1_mixer_below_valve.mat')
pos = [2 1 1 3 3 2];
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
c = parula(4);
clear l
for i = 1:width(data.MFC500)
	[y,x] = hist(data.MFC500(i,30e4:10:55e4),50);
	y = y/sum(y);
	plot(x*100,y,'Color',c(pos(i),:))
end
xlabel('Flow Rate (mL/min)')
subplot(1,2,2), hold on
c = parula(4);
for i = 1:width(data.PID)
	[y,x] = hist(data.PID(i,30e4:10:55e4),50);
	y = y/sum(y);
	l(pos(i)) = plot(x,y,'Color',c(pos(i),:));
end
xlabel('Stimulus (V)')
legend(l,{'Pos1','Pos2','Pos3'},'location','northeast');
PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% This is much better! The distribution doesn't dramatically change with the PID position, indicating that airstream is more well-mixed than previously. The errant distribution is probably due to some initial transient instead of some systemic error. 

%% Mean Shifted Gaussians
% We now attempt to generate a variety of Gaussians using this configuration: 

load('/local-data/DA-paper/MSG-construction/2015_07_22_MSG_mixer_below_valve.mat')

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(length(data)+1);
L = {};
l = [];
for i = 1:length(data)
	if ~isempty(data(i).PID)
		hist_this = data(i).PID(:,20e4:10:55e4)';
		xx =  linspace(min(min(hist_this)),max(max(hist_this)),50);
		y = NaN(width(hist_this),50);
		for j = 1:width(hist_this)
			y(j,:) = hist(hist_this(:,j),xx);
			y(j,:) = y(j,:)/sum(y(j,:));
		end
		l = [l errorShade(xx,mean(y),sem(y),'Color',c(i,:),'LineWidth',5)];
		L = [L ControlParadigm(i).Name]; 
	end
end
xlabel('Stimulus (V)')
set(gca,'XLim',[0 3.5])
legend(l,L)

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
