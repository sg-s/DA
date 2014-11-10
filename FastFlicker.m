% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Fast Flickering Odor Stimuli with different means
% In this project, we want to deliver a fast flickering odor stimulus that is approximately gaussian distributed. We want to vary the mean independently of the standard deviation. 

% internal housekeeping: determine if being called by publish or not
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end


%% Driving the MFCs as fast as possible
% Following a conversation with the very helpful Alicat staff, I changed the PID parameters of the 200mL/min MFC (P->500 and D->6000 from a default of P=100 and D=5000) to speed up the MFC dynamics. I then gave the MFC a control signal that took a new, uniformly distributed value every 10ms. 

%%
% The following figure shows what the control signal, and 3 repetitions of the MFC look like:


load('/local-data/DA-paper/fast-flicker/pid/2014_11_07_fast_flicker_single_MFC_2ac_10ms.mat')


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

time = 1e-4*(1:length(data(4).PID));

subplot(2,1,1), hold on
plot(time,ControlParadigm(4).Outputs(1,:),'k')
set(gca,'XLim',[30 35],'YLim',[-0.1 5])
ylabel('Control Signal (V)')

subplot(2,1,2), hold on
plot(time,data(4).MFC200')

set(gca,'XLim',[30 35],'YLim',[-0.1 5])
xlabel('Time (s)')
ylabel('Flow Signal (V)')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end



%%
% Apart from an inability to follow the control signal, it looks like the MFC is doing different things on different trials. The r-square between the first and the third trial is:

disp(rsquare(data(4).MFC200(1,200000:500000),data(4).MFC200(3,200000:500000)))

%%
% which shows that this is not good.


%% Driving the MFCs with a 50ms switching time
% Since the 10ms signal is too much for the MFC to handle, we drive it with a 50ms signal:

clear data
load('/local-data/DA-paper/fast-flicker/pid/2014_11_07_fast_flicker_single_MFC_2ac.mat')
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

time = 1e-4*(1:length(data(4).PID));

subplot(2,1,1), hold on
plot(time,ControlParadigm(4).Outputs(1,:),'k')
set(gca,'XLim',[30 35],'YLim',[-0.1 5])
ylabel('Control Signal (V)')

subplot(2,1,2), hold on
plot(time,data(4).MFC200')

set(gca,'XLim',[30 35],'YLim',[-0.1 5])
xlabel('Time (s)')
ylabel('Flow Signal (V)')


PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% Now the r-square between the 1st and 3rd trials is:
disp(rsquare(data(4).MFC200(1,200000:500000),data(4).MFC200(3,200000:500000)))

%%
% Now the r-square between the 2nd and 3rd trial is:
disp(rsquare(data(4).MFC200(2,200000:500000),data(4).MFC200(3,200000:500000)))

%%
% OK, what does the PID look like? The following figure shows all five cases: in each, the mean of the signal is increased compared to the last one. 

f=figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on

for i = 1:length(data)
	subplot(2,4,i); hold on
	plot(time,data(i).PID')
	set(gca,'XLim',[20 30],'YLim',[-0.1 7])
	title(ControlParadigm(i).Name)
end
xlabel('Time (s)')
ylabel('PID (V)')

PrettyFig;
if being_published
	snapnow
	delete(f)
end

%%
% What do the distributions look like?

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(length(data));
for i = 1:length(data)
	for j = 1:width(data(i).PID)
		[y,x] = hist(data(i).PID(j,200000:500000),50);
		plot(x,y,'Color',c(i,:))
	end
end

xlabel('PID (V)')
ylabel('Count')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% While the distributions are consistent from trial-to-trial, they in no way look Gaussian, nor do they shift along the X-axis in a nice manner. What is we drive the MFC with Gaussian instead of uniform noise? 




figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
c = parula(length(data));
for i = 1:length(data)
	for j = 1:width(data(i).PID)
		[y,x] = hist(data(i).PID(j,200000:500000),50);
		plot(x,y,'Color',c(i,:))
	end
end
title('MFC driven by uniform noise')
xlabel('PID (V)')
ylabel('Count')

subplot(1,2,2), hold on
clear data
load('/local-data/DA-paper/fast-flicker/pid/2014_11_07_gaussian_flicker_single_MFC_2ac_50ms.mat')
c = parula(length(data));
for i = 1:length(data)
	for j = 1:width(data(i).PID)
		[y,x] = hist(data(i).PID(j,200000:500000),50);
		plot(x,y,'Color',c(i,:))
	end
end

title('MFC driven by Gaussian noise')
xlabel('PID (V)')
ylabel('Count')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%% A static model to understand what is going on
% How do distributions of command voltages translate into distributions of flow signals to distributions of PID signals? In principle, we can work out the functional relationships between all these distributions. 



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
