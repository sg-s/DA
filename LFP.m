% LFP.m
% 
% created by Srinivas Gorur-Shandilya at 10:04 , 18 May 2015. Contact me at http://srinivas.gs/contact/
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

%% Measurements of the Local Field Potential 
% In this document, we attempt to measure the LFP. 

load('/local-data/DA-paper/LFP/2015_05_15_CS_F1_ab3_8_LFP_low_DC_gain1000.mat')

%% LFP response to pulses
% In the following figure, we present pulses of odor and record the LFP in response to these: 

figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
ax(1) = subplot(2,1,1); hold on
ax(2) = subplot(2,1,2); hold on
time = 1e-4*(1:length(data(10).PID));
plot(ax(1),time,mean2(data(10).PID),'b')
plot(ax(2),time,mean2(data(10).voltage),'b')
plot(ax(1),time,mean2(data(11).PID),'r')
plot(ax(2),time,mean2(data(11).voltage),'r')
ylabel(ax(2),'LFP (mV)')
ylabel(ax(1),'Odor stimulus (V)')
xlabel(ax(2),'Time (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Long-timescale changes in LFP in response to odor onsets
% When we start a 60-second odor flicker, we observe that the LFP changes very slowly in response to the overall increase in the odor stimulus: 

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,1,1), hold on
time = 1e-4*(1:length(data(19).PID));
plot(time,mean2(data(19).PID([1 5 10],:)),'k')
ylabel('PID (V)')
subplot(2,1,2), hold on

V= mean2(data(19).voltage([1 5 10],:));
plot(time,V,'k')
xlabel('Time (s)')
ylabel('LFP (V)')

[~,loc] = min(V);
ft=fittype('a*exp(-x./b1) +c');
fo = fitoptions(ft);
fo.Robust = 'on';
fo.StartPoint = [1.49 35 -.08];
ff = fit(time(1:end-loc+1)',-V(loc:end)',ft,fo);
l = plot(time(loc:end),-ff(time(1:end-loc+1)),'r');
legend(l,strcat('\tau=',oval(ff.b1),'s'),'Location','southeast')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

data(19).voltage(4,:) = []; % this trial is to be abandoned. 
data(19).PID(4,:) = [];

% Because there is a slow variation, and we're not really interested in it, we filter it out and plot the data. In the following figure, we remove spikes from the LFP, and then filter it with a high-pass filter with a 5 second cutoff (because we're not interested in any feature with a time scale longer than 5 seconds.) The LFP is shown together with the stimulus and the firing rates. 

c = parula(width(data(19).voltage)+1);

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(3,1,1), hold on
time = 1e-4*(1:length(data(19).PID));
for i = 1:width(data(19).PID)
	plot(time(1:20:end),data(19).PID(i,1:20:end),'Color',c(i,:))
end
ylabel('PID (V)')
set(gca,'XLim',[20 60])

time = time(1:10:end);
subplot(3,1,2), hold on
for i = 1:width(data(19).PID)
	this_LFP = data(19).voltage(i,:);
	% filter to remove spikes
	[~,this_LFP] = filter_trace(this_LFP);

	% now pass it through a high pass filter with a 5 second cutoff to remove slow fluctuations:
	this_LFP = this_LFP(1:10:end);
	temp = filtfilt(ones(5e3,1)/5e3,1,this_LFP);

	plot(time,this_LFP-temp,'Color',c(i,:))
end
ylabel('LFP (a.u.)')
set(gca,'XLim',[20 60],'YLim',[-.3 .3])


subplot(3,1,3), hold on
[f,t] = spiketimes2f(spikes(19).A,1e-4*(1:length(spikes(19).A)),1e-3);
f(:,4) = [];
for i = 1:width(f)
	plot(t,f(:,i),'Color',c(i,:))
end
ylabel('Firing Rate (Hz)')
set(gca,'XLim',[20 60],'YLim',[0 60])
xlabel('Time (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% Linear Models for LFP
% We now back out linear filters for the LFP on a trial-by-trial basis. 

K = zeros(601,width(data(19).PID));
for i = 1:width(data(19).PID)
	this_LFP = data(19).voltage(i,:);
	% filter
	% filter to remove spikes
	[~,this_LFP] = filter_trace(this_LFP);

	% now pass it through a high pass filter with a 5 second cutoff to remove slow fluctuations:
	this_LFP = this_LFP(1:10:end);
	temp = filtfilt(ones(5e3,1)/5e3,1,this_LFP);
	this_LFP = this_LFP - temp;
	this_LFP = this_LFP(20e3:end);
	time = 1e-3*(1:length(this_LFP));
	this_PID = data(19).PID(i,:);
	this_PID = this_PID(1:10:end);
	this_PID = this_PID(20e3:end);
	[K(:,i),~,filtertime] = FindBestFilter(this_PID,this_LFP,[],'regmax=1;','regmin=1;','filter_length=600;');
end

K = K(30:500,:);
filtertime = filtertime(30:500);

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[600 500]); hold on
plot(filtertime,K)
xlabel('Lag (ms)')
ylabel('Filter')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% How good are these predictions? The following figure shows an example trace with the best linear fit:


i = 2;
% filter
this_LFP = data(19).voltage(i,:);
% filter to remove spikes
[~,this_LFP] = filter_trace(this_LFP);

% now pass it through a high pass filter with a 5 second cutoff to remove slow fluctuations:
this_LFP = this_LFP(1:10:end);
temp = filtfilt(ones(5e3,1)/5e3,1,this_LFP);
this_LFP = this_LFP - temp;
this_LFP = this_LFP(20e3:end);
time = 1e-3*(1:length(this_LFP));
this_PID = data(19).PID(i,:);
this_PID = this_PID(1:10:end);
this_PID = this_PID(20e3:end);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(time+20,this_LFP,'k')
fp = convolve(time,this_PID,K(:,i),filtertime);
l = plot(time+20,fp,'r');
legend(strcat('r^2=',oval(rsquare(fp,this_LFP))))

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% How good are the predictions? The following figure shows the coefficient of determination for each of the trials:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

for i = 1:width(data(19).voltage)
	this_LFP = data(19).voltage(i,:);
	% filter to remove spikes
	[~,this_LFP] = filter_trace(this_LFP);

	% now pass it through a high pass filter with a 5 second cutoff to remove slow fluctuations:
	this_LFP = this_LFP(1:10:end);
	temp = filtfilt(ones(5e3,1)/5e3,1,this_LFP);
	this_LFP = this_LFP - temp;
	this_LFP = this_LFP(20e3:end);
	time = 1e-3*(1:length(this_LFP));
	this_PID = data(19).PID(i,:);
	this_PID = this_PID(1:10:end);
	this_PID = this_PID(20e3:end);
	fp = convolve(time,this_PID,K(:,i),filtertime);
	r2 = rsquare(fp,this_LFP);
	plot(i,r2,'k+')
end
xlabel('Trial')
ylabel('r^2')
set(gca,'YLim',[0 1])
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% 
% Maybe the prediction is poor because the LFP measurement is noisy. What if average all the LFPs, and then calculate the filters? 

mean_LFP = 0*this_LFP;
mean_PID = 0*this_PID;
for i = 1:width(data(19).voltage)
	this_LFP = data(19).voltage(i,:);
	% filter to remove spikes
	[~,this_LFP] = filter_trace(this_LFP);

	% now pass it through a high pass filter with a 5 second cutoff to remove slow fluctuations:
	this_LFP = this_LFP(1:10:end);
	temp = filtfilt(ones(5e3,1)/5e3,1,this_LFP);
	this_LFP = this_LFP - temp;
	this_LFP = this_LFP(20e3:end);
	time = 1e-3*(1:length(this_LFP));
	this_PID = data(19).PID(i,:);
	this_PID = this_PID(1:10:end);
	this_PID = this_PID(20e3:end);
	mean_PID = mean_PID + this_PID;
	mean_LFP = mean_LFP + this_LFP;
end
mean_LFP = mean_LFP/width(data(19).voltage);
mean_PID = mean_PID/width(data(19).voltage);
all_K = FindBestFilter(mean_PID,mean_LFP,[],'regmax=1;','regmin=1;','filter_length=600;');
all_K = all_K(30:500,:);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(time+20,mean_LFP,'k')
hold on
LFP_pred = convolve(time,mean_PID,all_K,filtertime);
l = plot(time+20,LFP_pred,'r');
legend(l,strcat('r^2=',oval(rsquare(LFP_pred,mean_LFP))));
xlabel('Time (s)')
ylabel('LFP (a.u.)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Linear Models for firing rate
% Why is the LFP so hard to predict with a linear model? What if we extract a filter for the firing rates directly from the PID?

K_Pf = zeros(1001,width(data(19).PID));
fp = NaN*f;
for i = 1:width(f)
	this_PID = data(19).PID(i,:);
	this_PID = this_PID(1:10:end);
	this_PID = this_PID(20e3:end);
	[K_Pf(:,i),~,filtertime] = FindBestFilter(this_PID,f(20e3:end,i),[],'regmax=1;','regmin=1;','filter_length=1000;','use_cache=0;');
end

K_Pf = K_Pf(70:570,:);
filtertime = filtertime(70:570);

for i = 1:width(f)
	this_PID = data(19).PID(i,:);
	this_PID = this_PID(1:10:end);
	this_PID = this_PID(20e3:end);
	fp(20e3:end,i) = convolve(time(20e3:end),this_PID,K_Pf(:,i),filtertime);
	temp = fit(fp(20e3:end-100,i),f(20e3:end-100,i),'poly1');
	fp(:,i) = temp(fp(:,i));	
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(t,mean2(f),'k')
l = plot(t,mean2(fp),'r');
legend(l,strcat('r^2=',oval(rsquare(mean2(f),mean2(fp)))))
set(gca,'XLim',[20 60])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% So, for some unknown reason, we can't predict the firing rates of the neuron very well. Something is weird about this, so it's hard to know what's going on with the LFP. 

%% Comparison of LFP and firing rates for odor and light
% How does the LFP get translated into firing? Here, we measure from flies expressing ReaChR in the ab3A neuron and activate that neuron with both light and odour. We then compare the LFP and the firing rate:

load('/local-data/DA-paper/reachr/2015_05_18_RR_F2_ab3_2_EA_2.mat')

haz_data =  [4 6 7 10];
figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[700 700]); hold on
for i = 1:length(haz_data)
	subplot(length(haz_data),2,2*(i-1)+1), hold on
	oktrials = setdiff(1:width(data(haz_data(i)).voltage),spikes(haz_data(i)).discard);
	V = [];
	for j = oktrials
		this_v= data(haz_data(i)).voltage(j,:);
		this_v(1:3e4) = [];
		[~,this_v] = filter_trace(this_v);
		this_v = this_v - mean(this_v(1:1e4));
		V = [this_v; V];
	end
	V = mean2(V);
	V = V(1:10:end);
	t = 1e-3*(1:length(V));
	plot(t,V)
	if i == 1
		title('LFP (mV)')
	end
	if any(strfind(ControlParadigm(haz_data(i)).Name,'Odour'))
		ylabel('Odour')
	else
		ylabel('Light')
	end

	subplot(length(haz_data),2,2*(i-1)+2), hold on
	[f,t]=(spiketimes2f(spikes(haz_data(i)).A));
	f = mean2(f);
	t= t-3;
	f(t<0)=[];
	t(t<0)=[];
	plot(t,f);
	set(gca,'YLim',[0 150])
	if i == 1
		title('Firing Rate (Hz)')
	end
end
xlabel('Time (s)')


PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% A trivial explanation for the discrepancy between firing rate/LFP for light and odour could be that the LFP is not a good measure of the total transmembrane current in the neuron. We could pick up only the current at the dendrite. Since we expect ReaChR to be localised everywhere on the neuron, and since it is fair to assume that the light penetrates the tissue very well, the LFP could be artifically small since we can only measure the transmembrance current in the dendrite, and not, say, in the cell body. 

%%
% A more interesting explanation would be that for the same transmembrane current, the neuron's firing machinery has different gains, based on whether the current is due to ORs or due to ReaChR. But this is improbable, and even if so, is very hard to show. 


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


