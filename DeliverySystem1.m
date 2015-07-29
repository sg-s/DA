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


%% Generalised NLN Model to build perfect distributions. 
% The broad goal of this section is to build a detailed, parametric model of the delivery system, and fit this to data. Then, using this model, we numerically optimise the distribution of control signals so that we get nice, Gaussian-distributed outputs from the delivery system. This effort is summarised in the following figure:
% 
% <</Users/sigbhu/code/da/images/2015-07-24-fig1.png>>
%
% 
% <</Users/sigbhu/code/da/images/2015-07-24-fig2.png>>
%


%%
% First, we fit our model from the control signals given to the MFC to the actual MFC flows.

fit_me = data(findData(data));
cp = ControlParadigm(findData(data));
for i = 1:length(fit_me)
	fit_me(i).stimulus = cp(i).Outputs(1,20e4:10:55e4);
	fit_me(i).response = mean2(fit_me(i).MFC500(:,20e4:10:55e4));
	fit_me(i).response(1:1e3) = NaN;
end

global p1
p1 = FitModel2Data(@pDeliverySystem,fit_me,'nsteps',0);

figure('outerposition',[0 0 1300 700],'PaperUnits','points','PaperSize',[1300 700]); hold on
for i = 1:length(fit_me)
	subplot(2,3,i), hold on
	t = 1e-3*(1:length(fit_me(i).response)) + 20;
	plot(t,fit_me(i).response,'k');
	pred = pDeliverySystem(fit_me(i).stimulus,p1);
	plot(t,pred,'r');
	set(gca,'XLim',[30 40],'YLim',[0.2 2])
	xlabel('Time (s)')
	ylabel('Stimulus (V)')
end

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end


% Now we fit the same generalised NLN model (with a different set of parameters) to account for the transformation from the MFC flows to the PID. 

fit_me = data(findData(data));
for i = 1:length(fit_me)
	fit_me(i).stimulus = mean2(fit_me(i).MFC500(:,20e4:10:55e4));
	fit_me(i).stimulus = (100*fit_me(i).stimulus)./(100*fit_me(i).stimulus+2000);
	fit_me(i).response = mean2(fit_me(i).PID(:,20e4:10:55e4));
	fit_me(i).response(1:1e3) = NaN;
end


global p2
p2 = FitModel2Data(@pDeliverySystem,fit_me,'nsteps',0);

figure('outerposition',[0 0 1300 700],'PaperUnits','points','PaperSize',[1300 700]); hold on
for i = 1:length(fit_me)
	subplot(2,3,i), hold on
	t = 1e-3*(1:length(fit_me(i).response)) + 20;
	plot(t,fit_me(i).response,'k');
	pred = pDeliverySystem(fit_me(i).stimulus,p2);
	plot(t,pred,'r');
	set(gca,'XLim',[30 40],'YLim',[.2 2.5])
	xlabel('Time (s)')
	ylabel('Stimulus (V)')
end

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% This model seems to do pretty well. How well does this model perform in estimating the distributions of the stimulus? 

figure('outerposition',[0 0 1300 700],'PaperUnits','points','PaperSize',[1300 700]); hold on
for i = 1:length(fit_me)
	subplot(2,3,i), hold on
	h = fit_me(i).response(1e3:end);
	[y,x] = hist(h,0:1e-2:5);
	y = y/sum(y);
	plot(x,y,'k')
	pred = pDeliverySystem(fit_me(i).stimulus,p2);
	[y,x] = hist(pred(1e3:end),0:1e-2:5);
	y = y/sum(y);
	plot(x,y,'r')
	xlabel('Stimulus (V)')
	set(gca,'XLim',[0.125 2.6])
end

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% While it's not perfect, it's pretty damn good. A better fit to the data would probably need a mechanistic understanding of how the delivery system works, and how odour moves from the liquid phase to the gas phase, which I don't have. 

%%
% We now specify the target stimulus distributions that we eventually want to produce with our delivery system:

clear p
p.   mu1= 0.5;
p.sigma1= 0.1;
p.   mu2= 0;
p.sigma2= 0;
p.  xmin= 0;
p.  xmax= 1;
x= 0:1e-2:5; 


c = parula(7);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:6
	ty = dist_gauss2(x,p);
	p.mu1 = p.mu1+.4;
	plot(x,ty,'Color',c(i,:))
end
xlabel('Stimulus (V)')

PrettyFig()

if being_published
	snapnow
	delete(gcf)
end

%%
% For each of these distributions, we use our detailed model of the delivery system to numerically optimise the distribution of control signals to get the best possible output distribution. The following figure shows the predicted output after optimising the input distributions: 


% the following code was actually used to fit 
% clear p
% p.   mu1= 0.5;
% p.sigma1= 0.1;
% p.   mu2= 0;
% p.sigma2= 0;
% p.  xmin= 0;
% p.  xmax= 1;
% x= 0:1e-2:5; 

% haz_data = [findData(data) 8];

% for i = 1:6
% 	clear d
% 	d.stimulus = x;
% 	d.response = hist(ControlParadigm(haz_data(i)).Outputs(1,20e4:10:55e4),x);
% 	% seed with fit to ansatz distribution
% 	pp = FitModel2Data(@dist_gamma2,d);

% 	clear d
% 	ty = dist_gauss2(x,p);
% 	d.response = ty/max(ty);
% 	d.stimulus = x;

% 	FitModel2Data(@BestDistribution,d,'UseParallel',false,'nsteps',300,'p0',pp);
% 	p.mu1 = p.mu1+.4;
% end


% recover from cache
clear p
p.   mu1= 0.5;
p.sigma1= 0.1;
p.   mu2= 0;
p.sigma2= 0;
p.  xmin= 0;
p.  xmax= 1;

clear pp
py = NaN(length(x),6);
for i = 1:6
	ty = dist_gauss2(x,p);
	d.response = ty/max(ty);
	d.stimulus = x;

	temp = FitModel2Data(@BestDistribution,d,'UseParallel',false,'nsteps',0);
	pp(i) = orderfields(temp);
	% also get the actual PID distributions
	py(:,i) = BestDistribution(0,pp(i));

	p.mu1 = p.mu1+.4;
end


clear p
p.   mu1= 0.5;
p.sigma1= 0.1;
p.   mu2= 0;
p.sigma2= 0;
p.  xmin= 0;
p.  xmax= 1;

c = parula(7);
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
for i = 1:6
	ty = dist_gauss2(x,p);
	p.mu1 = p.mu1+.4;
	plot(x,ty,'Color',c(i,:))
	plot(x,py(:,i),'Color',c(i,:))
end
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
