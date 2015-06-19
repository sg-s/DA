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


%% LFP responses to a flickering odour stimulus
% In the following figure we analyse the response of the LFP from two neurons to a flickering odour stimulus. The following figure shows the stimulus and the PID from the dataset. We bandpass the LFP to remove slow fluctuations that we don't care about, and remove spikes. 

load('/local-data/DA-paper/LFP/2015_06_16_RR_F4_ab3_11_EA.mat')
PID = data(22).PID;
LFP = data(22).voltage;
time = 1e-4*(1:length(LFP));
fA = spiketimes2f(spikes(22).A,time,1e-3);
t = 1e-3*(1:length(fA));
rm_this=find(sum(fA)==0);
PID(rm_this,:) = [];
LFP(rm_this,:) = [];
fA(:,rm_this) = [];
PID = PID(:,1:10:end);
LFP = LFP(:,1:10:end);
orn = ones(width(PID),1);

load('/local-data/DA-paper/LFP/2015_06_16_RR_F4_ab3_12_EA.mat')
PID2 = data(22).PID;
LFP2 = data(22).voltage;
time = 1e-4*(1:length(LFP2));
fA2 = spiketimes2f(spikes(22).A,time,1e-3);
PID2 = PID2(:,1:10:end);
LFP2 = LFP2(:,1:10:end);

rm_this=find(sum(fA2)==0);
PID2(rm_this,:) = [];
LFP2(rm_this,:) = [];
fA2(:,rm_this) = [];

PID = [PID; PID2];
LFP = [LFP; LFP2];
fA = [fA fA2];
orn = [orn; 2*ones(width(PID2),1)];

% bandpass
for i = 1:width(PID)
	LFP(i,:) = filter_trace(LFP(i,:),5000,10);
end

% remove mean 
for i = 1:width(PID)
	LFP(i,:) = LFP(i,:) -  mean(LFP(i,:));
end

time = 1e-3*(1:length(LFP));

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
l(1) = errorShade(time,mean2(PID(orn==1,:)),std(PID(orn==1,:)),'Color',[1 0 0],'SubSample',50);
l(2) = errorShade(time,mean2(PID(orn==2,:)),std(PID(orn==2,:)),'Color',[0 0 1],'SubSample',50);
set(gca,'XLim',[10 40])
ylabel('PID (V)')
legend(l,{'ORN 1','ORN 2'})
subplot(2,1,2), hold on
errorShade(time,mean2(LFP(orn==1,:)),std(LFP(orn==1,:)),'Color',[1 0 0]);
errorShade(time,mean2(LFP(orn==2,:)),std(LFP(orn==2,:)),'Color',[0 0 1]);
set(gca,'XLim',[10 40])
xlabel('Time (s)')
ylabel('LFP (norm)')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% This is a very clean dataset, and we see that the LFP is fairly reproducible, even between neurons. In this treatment, we ignore the absolute value of the LFP, by removing the mean and dividing through by the standard deviation. 

%% Linear filters
% Can a simple linear filter predict the LFP or the firing rate? In this section, we extract linear filters on a trial-wise basis for all combinations: from the PID to the LFP, from the PID to the firing rate, and from the LFP to the firing rate. 

% PID -> LFP
K_PL = zeros(width(PID),601);
offset = -100;
for i = 1:width(PID)
	K_PL(i,:) = FitFilter2Data(PID(i,1e4-offset:5e4-offset),LFP(i,1e4:5e4),[],'filter_length=600;','reg=1;');
end
K_PL = K_PL(:,50:550);
filtertime = 1e-3*(1:length(K_PL));
filtertime = filtertime  + (offset+50)*1e-3;

% LFP -> firing rate
K_Lf = zeros(width(PID),601);
offset = -100;
for i = 1:width(PID)
	K_Lf(i,:) = FitFilter2Data(LFP(i,1e4-offset:5e4-offset),fA(1e4:5e4,i),[],'filter_length=600;','reg=1;');
end
K_Lf = K_Lf(:,50:550);

% PID -> firing rate
K_Pf = zeros(width(PID),601);
offset = -100;
for i = 1:width(PID)
	K_Pf(i,:) = FitFilter2Data(PID(i,1e4-offset:5e4-offset),fA(1e4:5e4,i),[],'filter_length=600;','reg=1;');
end
K_Pf = K_Pf(:,50:550);


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
clear l
l(1) = errorShade(filtertime,mean2(K_PL(orn==1,:)),std(K_PL(orn==1,:)),'Color',[1 0 0]);
l(2) = errorShade(filtertime,mean2(K_PL(orn==2,:)),std(K_PL(orn==2,:)),'Color',[0 0 1]);
legend(l,{'ORN 1','ORN 2'},'Location','southeast')
xlabel('Filter Lag (s)')
ylabel('Filter Amplitude')
title('PID \rightarrow LFP')

subplot(1,3,2), hold on
l(1) = errorShade(filtertime,mean2(K_Lf(orn==1,:)),std(K_Lf(orn==1,:)),'Color',[1 0 0]);
l(2) = errorShade(filtertime,mean2(K_Lf(orn==2,:)),std(K_Lf(orn==2,:)),'Color',[0 0 1]);
legend(l,{'ORN 1','ORN 2'},'Location','southeast')
xlabel('Filter Lag (s)')
title('LFP \rightarrow Firing Rate')

subplot(1,3,3), hold on
l(1) = errorShade(filtertime,mean2(K_Pf(orn==1,:)),std(K_Pf(orn==1,:)),'Color',[1 0 0]);
l(2) = errorShade(filtertime,mean2(K_Pf(orn==2,:)),std(K_Pf(orn==2,:)),'Color',[0 0 1]);
legend(l,{'ORN 1','ORN 2'},'Location','northeast')
xlabel('Filter Lag (s)')
title('PID \rightarrow Firing Rate','interpreter','tex')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% How well do these filters perform? To assess this, we make linear predictions on a trial-wise basis and then compare the linear predictions to the data. 

LFP_pred = LFP*NaN;
for i = 1:length(orn)
	LFP_pred(i,:) = convolve(time,PID(i,:),K_PL(i,:),filtertime);
	r(i) = rsquare(LFP_pred(i,1e4:4e4),LFP(i,1e4:4e4));
end

figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
clear l
l(1) = errorShade(time,mean2(LFP)/std(mean2(LFP)),std(LFP),'Color',[0 0 0],'SubSample',10);
l(2) = errorShade(time,mean2(LFP_pred)/std(nonnans(mean2(LFP_pred))),std(LFP_pred),'Color',[1 0 0],'SubSample',10);
legend(l,{'mean LFP','Linear Prediction'})
xlabel('Time (s)')
ylabel('LFP (norm)')
set(gca,'XLim',[10 40])

subplot(2,1,2), hold on
plot(find(orn==1),r(find(orn==1)),'r-+')
plot(find(orn==2),r(find(orn==2)),'b-+')
legend({'ORN 1','ORN 2'},'Location','southeast')
set(gca,'YLim',[0 1])
xlabel('Trial #')
ylabel('r^2')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



%% Comparison of LFP and firing rates for odor and light
% How does the LFP get translated into firing? Here, we measure from flies expressing ReaChR in the ab3A neuron and activate that neuron with both light and odour. We then compare the LFP and the firing rate:

load('/local-data/DA-paper/LFP/2015_05_18_RR_F2_ab3_2_EA_2.mat')

haz_data =  [4 6 7 10];
figure('outerposition',[0 0 700 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
for i = 1:length(haz_data)
	subplot(length(haz_data),2,2*(i-1)+1), hold on
	oktrials = setdiff(1:width(data(haz_data(i)).voltage),spikes(haz_data(i)).discard);
	V = [];
	for j = oktrials
		this_v= data(haz_data(i)).voltage(j,:);
		this_v(1:3e4) = [];
		[~,this_v] = filter_trace(this_v,1000,Inf);
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


