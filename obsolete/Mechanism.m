% Mechanism.m
% in this document, we attempt to understand the mechanism behind gain adaptation
% 
% created by Srinivas Gorur-Shandilya at 1:15 , 06 April 2015. Contact me at http://srinivas.gs/contact/
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


%% What is the mechanism behind gain adaptation?
% Two (non-exclusive) possibilities are that the gain is controlled at the receptor level, and that the gain is controlled at the firing machinery. To disambiguate the two, we record the neuron's responses to a flickering odour stimulus. We then repeat the experiment, but also drive the neuron optically through light activating ReaChR channels. 

%% The Microscope light is effective in activating the neuron
% In fact, it is far more effective than the specially built LED light we built for the purpose. The following figure shows the set up illuminated by the high-powered LEDs:
% 
% <</Users/sigbhu/code/da/images/amber-led.jpg>>
%

%%
% Despite being really bright, the power at 591nm at the location of the fly from these LEDs isn't even close to what the microscope light can crank out:

LightLevels = [287 315 315; 297 323 308; 294 313 316; 304 309 313; 246 241 241 ; 127 162 165; 771 727 725; 1280 1320 1250; 2210 2060 2070];

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
errorbar(mean(LightLevels'),std(LightLevels'),'kx')
set(gca,'XTickLabel',{'LED:5V','LED:4V','LED:3V','LED:2V','LED:1V','M: 2V','M: 3V','M: 3.5V','M: 4V'})
set(gca,'XTick',[1:9],'XTickLabelRotation',45)
ylabel('Power @ 591nm (\muW)')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% The following figure shows a trace where the microscope light was turned on around the 10s mark. 

load('/local-data/DA-paper/reachr/2015_04_03_ReaChR_F4_ab3_1_EA.mat')
figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
time = 1e-4*(1:length(data(6).PID));
plot(time,data(6).voltage,'k')
plot(time(logical(spikes(6).A)),data(6).voltage(logical(spikes(6).A)),'ro');
plot(time(logical(spikes(6).B)),data(6).voltage(logical(spikes(6).B)),'bo');
set(gca,'YLim',[-.1 .1])
xlabel('Time (s)')
ylabel('Voltage (\muV)')

PrettyFig('plw=1;');

if being_published
	snapnow
	delete(gcf)
end

%% ORN responses with and without light stimulation
% In the following figure, we show data from 1 sensillum where we record responses to a flickering odour stimulus with and without light stimulus. 

load('/local-data/DA-paper/reachr/2015_04_03_ReaChR_F3_ab3_3_EA.mat')

fA_no_light = spiketimes2f(spikes(2).A(2:3,:),time,1e-3);
fA = spiketimes2f(spikes(8).A,time,1e-3);
tA = 1e-3*(1:length(fA));
PID = fA;
for i = 1:width(PID)
	PID(:,i) = interp1(time,data(8).PID(i,:),tA);
end
PID_no_light = fA_no_light;
for i = 1:width(PID_no_light)
	PID_no_light(:,i) = interp1(time,data(2).PID(i+1,:),tA);
end

%%
% First, we show the distributions of the stimulus and the response without light (blue) vs with light (red). 
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ax(1) = subplot(1,2,1); hold on
xlabel('Odour Stimulus (V)')
ax(2) = subplot(1,2,2); hold on
xlabel('ORN response (Hz)')
for i = 1:width(PID)
	[hy,hx] = hist(PID(:,i),50);
	plot(ax(1),hx,hy,'r')
	[hy,hx] = hist(fA(:,i),50);
	plot(ax(2),hx,hy,'r')
end
for i = 1:width(PID_no_light)
	[hy,hx] = hist(PID_no_light(:,i),50);
	plot(ax(1),hx,hy,'b')
	[hy,hx] = hist(fA_no_light(:,i),50);
	plot(ax(2),hx,hy,'b')
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% Now we extract filters for the odour and the no-odour case and compare them, and also compare the response vs. the linear prediction. 

K = [];
fp = fA;
for i = 1:width(PID)
	[this_K, ~, filtertime_full] = FindBestFilter(PID(:,i),fA(:,i),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*mean(diff(tA));
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	this_K = this_K/max(this_K);
	K = [K;this_K];
	fp(:,i) = convolve(tA,PID(:,i),this_K,filtertime);
end

K_no_light = [];
fp_no_light = fA_no_light;
for i = 1:width(PID_no_light)
	[this_K, ~, filtertime_full] = FindBestFilter(PID_no_light(:,i),fA_no_light(:,i),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*mean(diff(tA));
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	this_K = this_K/max(this_K);
	K_no_light = [K_no_light;this_K];
	fp_no_light(:,i) = convolve(tA,PID_no_light(:,i),this_K,filtertime);
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
shadedErrorBar(filtertime,mean2(K),std(K),{'Color','r'})
shadedErrorBar(filtertime,mean2(K_no_light),std(K_no_light),{'Color','b'})
xlabel('Lag (s)')
ylabel('Filter amplitude (norm)')

subplot(1,2,2), hold on
x = mean2(fp_no_light);
y = mean2(fA_no_light);
ss = 50;
plot(x(1:ss:end),y(1:ss:end),'b.')
x = mean2(fp);
y = mean2(fA);
plot(x(1:ss:end),y(1:ss:end),'r.')
xlabel('K \otimes s')
ylabel('Response (Hz)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Another neuron
% We now look at data from another neuron, where we do the same experiment. Here we have three cases: no light, and 3, 3.5 of microscope light. (we have one trial of 4V, but the neuron starts pinching and drifts)


load('/local-data/DA-paper/reachr/2015_04_03_ReaChR_F1_ab3_1_EA.mat')
% orn response and PID
fA_no_light = spiketimes2f(spikes(2).A([1 3 4],:),time,1e-3);
PID_no_light = fA_no_light;
these=[1 3 4];
for i = 1:width(PID_no_light)
	PID_no_light(:,i) = interp1(time,data(2).PID(these(i),:),tA);
end

fA_3V = spiketimes2f(spikes(3).A(4:5,:),time,1e-3);
PID_3V = fA_3V;
for i = 1:width(PID_3V)
	PID_3V(:,i) = interp1(time,data(3).PID(3+i,:),tA);
end

fA_3_5V = spiketimes2f(spikes(3).A(6:9,:),time,1e-3);
PID_3_5V = fA_3_5V;
for i = 1:width(PID_3_5V)
	PID_3_5V(:,i) = interp1(time,data(3).PID(5+i,:),tA);
end

tA = 1e-3*(1:length(fA_3V));
c = parula(4);


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ax(1) = subplot(1,2,1); hold on
xlabel('Odour Stimulus (V)')
ax(2) = subplot(1,2,2); hold on
xlabel('ORN response (Hz)')
for i = 1:width(PID_no_light)
	[hy,hx] = hist(PID_no_light(:,i),50);
	plot(ax(1),hx,hy,'Color',c(1,:))
	[hy,hx] = hist(fA_no_light(:,i),50);
	plot(ax(2),hx,hy,'Color',c(1,:))
end
for i = 1:width(PID_3V)
	[hy,hx] = hist(PID_3V(:,i),50);
	plot(ax(1),hx,hy,'Color',c(2,:))
	[hy,hx] = hist(fA_3V(:,i),50);
	plot(ax(2),hx,hy,'Color',c(2,:))
end
for i = 1:width(PID_3_5V)
	[hy,hx] = hist(PID_3_5V(:,i),50);
	plot(ax(1),hx,hy,'Color',c(3,:))
	[hy,hx] = hist(fA_3_5V(:,i),50);
	plot(ax(2),hx,hy,'Color',c(3,:))
end

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

% We now compare the filters and the linear predictions:

K_no_light = [];
fp_no_light = fA_no_light;
for i = 1:width(PID_no_light)
	[this_K, ~, filtertime_full] = FindBestFilter(PID_no_light(:,i),fA_no_light(:,i),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*mean(diff(tA));
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	this_K = this_K/max(this_K);
	K_no_light = [K_no_light;this_K];
	fp_no_light(:,i) = convolve(tA,PID_no_light(:,i),this_K,filtertime);
end

K_3V = [];
fp_3V = fA_3V;
for i = 1:width(PID_3V)
	[this_K, ~, filtertime_full] = FindBestFilter(PID_3V(:,i),fA_3V(:,i),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*mean(diff(tA));
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	this_K = this_K/max(this_K);
	K_3V = [K_3V;this_K];
	fp_3V(:,i) = convolve(tA,PID_3V(:,i),this_K,filtertime);
end

K_3_5V = [];
fp_3_5V = fA_3_5V;
for i = 1:width(PID_3_5V)
	[this_K, ~, filtertime_full] = FindBestFilter(PID_3_5V(:,i),fA_3_5V(:,i),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
	filtertime_full = filtertime_full*mean(diff(tA));
	filtertime = 1e-3*(-200:900);
	this_K = interp1(filtertime_full,this_K,filtertime);
	this_K = this_K/max(this_K);
	K_3_5V = [K_3_5V;this_K];
	fp_3_5V(:,i) = convolve(tA,PID_3_5V(:,i),this_K,filtertime);
end




figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
shadedErrorBar(filtertime,mean2(K_no_light),std(K_no_light),{'Color',c(1,:)})
shadedErrorBar(filtertime,mean2(K_3V),std(K_3V),{'Color',c(2,:)})
shadedErrorBar(filtertime,mean2(K_3_5V),std(K_3_5V),{'Color',c(3,:)})

xlabel('Lag (s)')
ylabel('Filter amplitude (norm)')

subplot(1,2,2), hold on
x = mean2(fp_no_light);
y = mean2(fA_no_light);
ss = 50;
l(1)=plot(x(1:ss:end),y(1:ss:end),'.','Color',c(1,:));

x = mean2(fp_3V);
y = mean2(fA_3V);
l(2)=plot(x(1:ss:end),y(1:ss:end),'.','Color',c(2,:));

x = mean2(fp_3_5V);
y = mean2(fA_3_5V);
l(3)=plot(x(1:ss:end),y(1:ss:end),'.','Color',c(3,:));

legend(l,{'No light','3V','3.5V'})
xlabel('K \otimes s')
ylabel('Response (Hz)')

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

