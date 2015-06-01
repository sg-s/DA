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

% So we ignore for the time being trials where we start from no odour -> odour. 

%% LFP changes with flickering odor stimulus
% Here, we present a flickering odor stimulus, wait for some time while the slow changes in the LFP go away, and then record the LFP responses. In the following figure, we remove spikes from the LFP, mean subtract, and divide by the standard deviation:

do_these = [2 3 6 7 8 9 11];
c = parula(length(do_these)+1);

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,1,1), hold on
time = 1e-4*(1:length(data(19).PID));
for i = 1:length(do_these)
	plot(time,data(19).PID(do_these(i),:),'Color',c(i,:))
end
ylabel('PID (V)')
set(gca,'XLim',[20 60])

subplot(2,1,2), hold on
time = 1e-4*(1:length(data(19).PID));
for i = 1:length(do_these)
	this_LFP = data(19).voltage(do_these(i),:);
	% filter
	[~,this_LFP] = filter_trace(this_LFP);
	this_LFP = this_LFP - mean(this_LFP(2e5:end));
	this_LFP = this_LFP/std(this_LFP(2e5:end));
	plot(time,this_LFP,'Color',c(i,:))
end
ylabel('LFP (mV)')
set(gca,'XLim',[20 60])
xlabel('Time (s)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%% Linear Models for LFP
% We now back out linear filters for the LFP on a trial-by-trial basis. 

K = zeros(601,length(do_these));
for i = 1:length(do_these)
	this_LFP = data(19).voltage(do_these(i),:);
	% filter
	[~,this_LFP] = filter_trace(this_LFP);
	this_LFP = this_LFP - mean(this_LFP(2e5:end));
	this_LFP = this_LFP/std(this_LFP(2e5:end));
	time = 1e-4*(1:length(this_LFP));
	t = 20+  (1e-3:1e-3:40);
	this_LFP = interp1(time,this_LFP,t);
	this_PID = data(19).PID(do_these(i),:);
	this_PID = interp1(time,this_PID,t);
	this_PID = this_PID - mean(this_PID);
	this_PID = this_PID/std(this_PID);
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
this_LFP = data(19).voltage(do_these(i),:);
% filter
[~,this_LFP] = filter_trace(this_LFP);
this_LFP = this_LFP - mean(this_LFP(2e5:end));
this_LFP = this_LFP/std(this_LFP(2e5:end));
time = 1e-4*(1:length(this_LFP));
t = 20+  (1e-3:1e-3:40);
this_LFP = interp1(time,this_LFP,t);
this_PID = data(19).PID(do_these(i),:);
this_PID = interp1(time,this_PID,t);
this_PID = this_PID - mean(this_PID);
this_PID = this_PID/std(this_PID);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(t,this_LFP,'k')
fp = convolve(t,this_PID,K(:,end),filtertime);
plot(t,fp,'r');


PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% How good are the predictions? The following figure shows the coefficient of determination for each of the trials:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

for i = 1:length(do_these)
	this_LFP = data(19).voltage(do_these(i),:);
	% filter
	[~,this_LFP] = filter_trace(this_LFP);
	this_LFP = this_LFP - mean(this_LFP(2e5:end));
	this_LFP = this_LFP/std(this_LFP(2e5:end));
	time = 1e-4*(1:length(this_LFP));
	t = 20+  (1e-3:1e-3:40);
	this_LFP = interp1(time,this_LFP,t);
	this_PID = data(19).PID(do_these(i),:);
	this_PID = interp1(time,this_PID,t);
	this_PID = this_PID - mean(this_PID);
	this_PID = this_PID/std(this_PID);
	fp = convolve(t,this_PID,K(:,end),filtertime);
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
for i = 1:length(do_these)
	this_LFP = data(19).voltage(do_these(i),:);
	% filter
	[~,this_LFP] = filter_trace(this_LFP);
	this_LFP = this_LFP - mean(this_LFP(2e5:end));
	this_LFP = this_LFP/std(this_LFP(2e5:end));
	time = 1e-4*(1:length(this_LFP));
	t = 20+  (1e-3:1e-3:40);
	this_LFP = interp1(time,this_LFP,t);
	mean_LFP = mean_LFP + this_LFP;

	this_PID = data(19).PID(do_these(i),:);
	this_PID = interp1(time,this_PID,t);
	this_PID = this_PID - mean(this_PID);
	this_PID = this_PID/std(this_PID);
	mean_PID = mean_PID + this_PID;
end
mean_LFP = mean_LFP/length(do_these);
mean_PID = mean_PID/length(do_these);
all_K = FindBestFilter(this_PID,this_LFP,[],'regmax=1;','regmin=1;','filter_length=600;');
all_K = all_K(30:500,:);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(t,mean_LFP,'k')
hold on
fp = convolve(t,mean_PID,all_K,filtertime);
l = plot(t,fp,'r');
legend(l,strcat('r^2=',oval(rsquare(fp,mean_LFP))));
xlabel('Time (s)')
ylabel('LFP (norm)')
PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


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


