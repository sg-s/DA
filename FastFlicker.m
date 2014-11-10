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


%% A dynamic model to understand what's going on
% In this section, we build a detailed phenomenological model of the odor delivery system, to try to understand what's going on. We model all dynamical processes using LN models that are fit to the data, and assume that the PID measures instantaneous odor concentration that scales with the fractional flow through the odor. 
% 
% <</code/da/images/model.png>>
%

clear data
load('/local-data/DA-paper/fast-flicker/pid/2014_11_07_fast_flicker_single_MFC_2ac.mat')

% extract all MFC filters
if ~exist('K1','var')
	K1 = zeros(length(data)*width(data(1).PID),500);
	c= 1;
	for i = 1:length(data)
		for j = 1:width(data(i).PID)
			K1(c,:) = FitFilter2Data(ControlParadigm(i).Outputs(1,200000:10:500000),data(i).MFC200(j,200000:10:500000),[],'filter_length=499;');
			c = c+1;
			
		end
	end
end


% extract all non-linearities for MFC
if ~exist('f1','var')
	f1 = {};
	c= 1;
	for i = 1:length(data)
		for j = 1:width(data(i).PID)
			xdata = ControlParadigm(i).Outputs(1,200000:10:500000);
			xdata = filter(K1(c,:),1,xdata);
			ydata = data(i).MFC200(j,200000:10:500000);

			% crop it to lose NaNs
			ydata(isnan(xdata)) = [];
			xdata(isnan(xdata)) = [];

			xdata = xdata(:);
			ydata = ydata(:);

			f1{c} = fit(xdata,ydata,'poly6');
			c= c+1;
		end
	end
end

%%
% The following figure shows the filters extracted for the MFC response from the command signal, for each trial, colour coded by the control paradigm (colours match previous figures). The panel on the right shows the best-fit non linear functions (which are essentially straight lines). And the two panels below show the data (black) vs the prediction (red) for the lowest flow case and the highest flow case. 

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(3,2,1), hold on
c = parula(length(data));
filtertime = 1e-3*(1:length(K1));
for i = 1:width(K1)
	plot(filtertime,K1(i,:),'Color',c(ceil(i/3),:))
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (a.u.)')
title('MFC Filters')

subplot(3,2,2), hold on
c = parula(length(data));
filtertime = 1e-3*(1:length(K1));
cc= 1;
for i = 1:length(data)
	for j = 1:width(data(1).PID)
		xdata = ControlParadigm(i).Outputs(1,200000:10:500000);
		xdata = filter(K1(cc,:),1,xdata);
		xdata = sort(xdata);
		plot(xdata,f1{cc}(xdata),'Color',c(i,:));
		cc = cc+1;
	end
end
xlabel('Filter output (a.u.)')
ylabel('MFC Flow Rate (V)')
title('MFC Nonlinearities')


subplot(3,2,3:4), hold on
time = 1e-4*(1:length(data(1).PID));
plot(time,data(1).MFC200(3,:),'k')
xdata = ControlParadigm(1).Outputs(1,200000:10:500000);
xdata = filter(K1(3,:),1,xdata);
fp = f1{3}(xdata);
plot(time(200000:10:500000),fp,'r')
set(gca,'XLim',[20 30])
xlabel('Time (s)')
ylabel('Flow Signal (V)')

subplot(3,2,5:6), hold on
plot(time,data(7).MFC200(3,:),'k')
xdata = ControlParadigm(7).Outputs(1,200000:10:500000);
xdata = filter(K1(21,:),1,xdata);
fp = f1{21}(xdata);
plot(time(200000:10:500000),fp,'r')
set(gca,'XLim',[20 30])
xlabel('Time (s)')
ylabel('Flow Signal (V)')



PrettyFig;
if being_published
	snapnow
	delete(gcf)
end





% extract all PID filters
if ~exist('K2','var')
	K2 = zeros(length(data)*width(data(1).PID),500);
	c= 1;
	for i = 1:length(data)

		for j = 1:width(data(i).PID)
			x = data(i).MFC200(j,200000:10:500000);
			% convert into a flow
			x = x*40; % mL/min
			x = x./(x+2000);
			K2(c,:) = FitFilter2Data(x,data(i).PID(j,200000:10:500000),[],'filter_length=499;');
			c = c+1;

		end
	end
end

% extract all non-linearities for PID
if ~exist('f2','var')
	f2 = {};
	c= 1;
	for i = 1:length(data)
		for j = 1:width(data(i).PID)
			xdata = data(i).MFC200(j,200000:10:500000);
			% convert into a flow
			x = x*40; % mL/min
			x = x./(x+2000);
			xdata = filter(K2(c,:),1,xdata);
			ydata = data(i).PID(j,200000:10:500000);

			% crop it to lose NaNs
			ydata(isnan(xdata)) = [];
			xdata(isnan(xdata)) = [];

			xdata = xdata(:);
			ydata = ydata(:);

			f2{c} = fit(xdata,ydata,'poly6');
			c= c+1;
		end
	end
end


figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(3,2,1), hold on
c = parula(length(data));
filtertime = 1e-3*(1:length(K2));
for i = 1:width(K2)
	plot(filtertime,K2(i,:),'Color',c(ceil(i/3),:))
end
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude (a.u.)')
title('PID Filters')

subplot(3,2,2), hold on
c = parula(length(data));
filtertime = 1e-3*(1:length(K1));
cc= 1;
for i = 1:length(data)
	for j = 1:width(data(1).PID)
		xdata = data(i).MFC200(j,200000:10:500000);
		xdata = filter(K2(cc,:),1,xdata);
		xdata = sort(xdata);
		plot(xdata,f2{cc}(xdata),'Color',c(i,:));
		cc = cc+1;
	end
end
xlabel('Filter output (a.u.)')
ylabel('PID (V)')
title('PID Nonlinearities')


subplot(3,2,3:4), hold on
plot(time,data(1).PID(3,:),'k')
xdata = data(1).MFC200(3,200000:10:500000);
xdata = filter(K2(3,:),1,xdata);
fp = f2{3}(xdata);
plot(time(200000:10:500000),fp,'r')
set(gca,'XLim',[20 30])
xlabel('Time (s)')
ylabel('PID (V)')

subplot(3,2,5:6), hold on
plot(time,data(7).PID(3,:),'k')
xdata = data(7).MFC200(3,200000:10:500000);
xdata = filter(K2(21,:),1,xdata);
fp = f2{21}(xdata);
plot(time(200000:10:500000),fp,'r')
set(gca,'XLim',[20 30])
xlabel('Time (s)')
ylabel('PID (V)')



PrettyFig;
if being_published
	snapnow
	delete(gcf)
end


%% Wrapped Noise
% In this section, we choose our control signals by wrapping uniformly distributed noise around the extrema, and see what sort of distributions we get out of this. 

load('/local-data/DA-paper/fast-flicker/pid/2014_11_07_fast_flicker_single_MFC_2ac_50ms_wrapped_noise.mat')
figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,2,1), hold on
c = jet(length(data));
for i = 1:length(data)
	for j = 1:width(data(i).PID)
		[y,x] = hist(data(i).MFC200(j,200000:500000),50);
		plot(x,y,'Color',c(i,:))
	end
end

xlabel('MFC Flow (V)')
ylabel('Count')


subplot(2,2,2), hold on
c = jet(length(data));
for i = 1:length(data)
	for j = 1:width(data(i).PID)
		[y,x] = hist(data(i).PID(j,200000:500000),50);
		plot(x,y,'Color',c(i,:))
	end
end

xlabel('PID (V)')
ylabel('Count')

subplot(2,2,3:4), hold on
c = jet(length(data));
for i = 1:length(data)
	if ~isempty(data(i).PID)
		time = 1e-4*(1:length(data(i).PID));
		plot(time,mean(data(i).PID),'Color',c(i,:))
	end
end

set(gca,'XLim',[20 30])

xlabel('Time (s)')
ylabel('PID (V)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end



%% Dilution Noise
% In this section, we choose our control signals by sampling from a uniform distribution in the effective dilution (the flow through odor/total flow), and then back-calculate the command signal needed to drive the MFC. 

load('/local-data/DA-paper/fast-flicker/pid/2014_11_07_fast_flicker_single_MFC_2ac_50ms_dilution_noise.mat')
figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,2,1), hold on
c = jet(length(data));
for i = 1:length(data)
	for j = 1:width(data(i).PID)
		[y,x] = hist(data(i).MFC200(j,200000:500000),50);
		plot(x,y,'Color',c(i,:))
	end
end

xlabel('MFC Flow (V)')
ylabel('Count')


subplot(2,2,2), hold on
c = jet(length(data));
for i = 1:length(data)
	for j = 1:width(data(i).PID)
		[y,x] = hist(data(i).PID(j,200000:500000),50);
		plot(x,y,'Color',c(i,:))
	end
end

xlabel('PID (V)')
ylabel('Count')

subplot(2,2,3:4), hold on
c = jet(length(data));
for i = 1:length(data)
	if ~isempty(data(i).PID)
		time = 1e-4*(1:length(data(i).PID));
		plot(time,mean(data(i).PID),'Color',c(i,:))
	end
end

set(gca,'XLim',[20 30])

xlabel('Time (s)')
ylabel('PID (V)')


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
