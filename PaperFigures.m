% Paper Figures
% makes all the figures for the paper
% 
% created by Srinivas Gorur-Shandilya at 12:57 , 21 January 2015. Contact me at http://srinivas.gs/contact/
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

%% Figure 1

clearvars -except being_published
load('/local-data/DA-paper/carlotta/fig3/abc.mat')
figure('outerposition',[0 0 1100 700],'PaperUnits','points','PaperSize',[1100 700]); hold on
axes_handles = NaN(6,1);
for i = 1:6
	axes_handles(i) = subplot(2,3,i);
	hold on
end

% combine all spikes
all_spikes = [];
paradigm = [];
for i = 1:length(spikes)
	if isempty(spikes(i).discard)
		all_spikes = vertcat(all_spikes,spikes(i).A);
		paradigm = [paradigm i*ones(1,width(spikes(i).A))];
	else
		temp = spikes(i).A;
		rm_this = find(spikes(i).discard);
		rm_this(rm_this>width(spikes(i).discard)) = [];
		temp(rm_this,:) = [];
		all_spikes = vertcat(all_spikes,temp);
		paradigm = [paradigm i*ones(1,width(temp))];
	end
end

time = 1e-4*(1:length(all_spikes));

hash = DataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end

t = 1e-3*(1:length(fA));


% remove bizzare outliers
rm_this(1) = 8+ find(paradigm==6,1,'first');
rm_this(2) = 2+ find(paradigm==8,1,'first');
rm_this = [rm_this find(max(fA)==0)];
fA(:,rm_this) = [];
paradigm(rm_this) = [];

% combine all PID traces. the data is very spotty, so we have to be careful. 
all_PID = [];
paradigm_PID = [];
for i = 1:length(spikes)
	temp = data(i).PID;
	temp = temp(find(temp(:,1)),:);
	all_PID = vertcat(all_PID,temp);
	paradigm_PID = [paradigm_PID i*ones(1,width(temp))];
end
all_PID = all_PID';

% filter and subsample
for i = 1:width(all_PID)
	all_PID(:,i) = filter(ones(30,1)/30,1,all_PID(:,i));
end
all_PID = all_PID(1:10:end,:);

% remove bizzare outliers
rm_this = [];
rm_this(1) = find(paradigm_PID==9,1,'first');
rm_this(2) = find(paradigm_PID==10,1,'first');
all_PID(:,rm_this) = [];
paradigm_PID(rm_this) = [];

% plot PID
for i = 1:max(paradigm)
	plot_these = find(paradigm_PID==i);
	plot(axes_handles(1),t,mean2(all_PID(:,plot_these)));
end

set(axes_handles(1),'XLim',[.5 3],'YScale','log','YLim',[7e-2 15])
xlabel(axes_handles(1),'Time (s)')
ylabel(axes_handles(1),'Stimulus (V)')


% plot neuron responses
for i = 1:max(paradigm)
	plot_these = find(paradigm==i);
	plot(axes_handles(2),t,mean2(fA(:,plot_these)));
end

set(axes_handles(2),'XLim',[.5 3])
xlabel(axes_handles(2),'Time (s)')
ylabel(axes_handles(2),'Firing rate (Hz)')

% compute some stuff
nominal_stimulus_start = 1.05;
nominal_stimulus_stop = 1.6;
ttp.data = NaN(width(fA),1);
ttp.stimulus = NaN(width(fA),1);
peak.data = NaN(width(fA),1);
dga.data = NaN(width(fA),1);
for i = 1:width(fA)
	[m,loc]=max(fA(nominal_stimulus_start*1e3:nominal_stimulus_stop*1e3,i));
	ttp.data(i) = loc*1e-3;
	ttp.stimulus(i) = max(max(all_PID(:,paradigm_PID == paradigm(i))));
	peak.data(i) = m;
	dga.data(i) = mean(fA(1500:1600,i))/m;
end

% plot time to peak for the data
plot(axes_handles(4),ttp.stimulus,ttp.data,'k+')
set(axes_handles(4),'XScale','log','XLim',[5e-2 20])
xlabel(axes_handles(4),'Mean Stimulus (V)')
ylabel(axes_handles(4),'Time to peak (s)')

% plot dose-response for data
plot(axes_handles(5),ttp.stimulus,peak.data,'k+')
xlabel(axes_handles(5),'Mean Stimulus (V)')
set(axes_handles(5),'XScale','log','XLim',[5e-2 20])
ylabel(axes_handles(5),'Peak response (Hz)')

% plot degree of adaptation for data
plot(axes_handles(6),ttp.stimulus,dga.data,'k+')
xlabel(axes_handles(6),'Mean Stimulus (V)')
set(axes_handles(6),'XScale','log','XLim',[5e-2 20],'YLim',[0 1])
ylabel(axes_handles(6),'Degree of Adaptation')

% use LN model to make predictions
clear p
p.  tau1= 5.6925;
p.   K_n= 6.9533;
p.  tau2= 20.1989;
p.   K_A= 0.5832;
p.offset= 0.9780;
p. scale= 9.4446;
p. A= 230.9741;
p. k= 1.1985;
p. n= 0.7189;


K = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,1:1000);
K = K/max(K);

fp = NaN(length(fA),length(unique(paradigm)));
for i = 2:max(paradigm)
	this_pid = mean2(all_PID(:,paradigm_PID==i));
	this_pid = this_pid - mean(this_pid(5:1e3));
	this_pid(this_pid<0) =0;
	fp(:,i) = filter(K,1,this_pid);
	fp(:,i) = fp(:,i)/max(fp(:,i));
	fp(:,i) = fp(:,i)*hill(p,mean(ttp.stimulus(paradigm==i)));
end
fp(fp<0) = 0;
plot(axes_handles(3),t,fp,'r')
xlabel(axes_handles(3),'Time (s)')
set(axes_handles(3),'XLim',[.5 3])

plot(axes_handles(5),[5e-2:.1:20],hill(p,[5e-2:.1:20]),'r')

% compute metrics for LN model
[m,loc]=max(fp);
ttp.fp = loc*1e-3- nominal_stimulus_start;
dga.fp = mean(fp(1500:1600,:))./m;
dga.fp(1) = [];
x = sort(unique(ttp.stimulus));
x(1) = []; ttp.fp(1) = [];

% plot metrics for LN model
plot(axes_handles(4),x,ttp.fp,'r')
plot(axes_handles(6),x,dga.fp,'r')

PrettyFig('plw=1.5;','lw=1.5;','fs=14;')

if being_published
	snapnow
	delete(gcf)
end



%% Figure 2

%% Figure 3

%% Figure 4. 


clearvars -except being_published
fig_handle=figure('Units','pixels','outerposition',[82 5 971 851],'PaperUnits','points','PaperSize',[971 851],'Color','w','Toolbar','none');
clf(fig_handle); clear axes_handles
axes_handles(1)=axes('Units','pixels','Position',[63.825 723.35 638.25 85.1]);
axes_handles(2)=axes('Units','pixels','Position',[63.825 595.7 638.25 106.375]);
axes_handles(3)=axes('Units','pixels','Position',[63.825 510.6 638.25 63.825]);
axes_handles(4)=axes('Units','pixels','Position',[765.9 680.8 127.65 127.65]);
axes_handles(5)=axes('Units','pixels','Position',[765.9 510.6 127.65 127.65]);
axes_handles(6)=axes('Units','pixels','Position',[85.1 255.3 191.475 191.475]);
axes_handles(7)=axes('Units','pixels','Position',[85.1 42.55 191.475 170.2]);
axes_handles(8)=axes('Units','pixels','Position',[382.95 255.3 191.475 191.475]);
axes_handles(9)=axes('Units','pixels','Position',[382.95 42.55 191.475 170.2]);
axes_handles(10)=axes('Units','pixels','Position',[702.075 255.3 191.475 191.475]);
axes_handles(11)=axes('Units','pixels','Position',[702.075 42.55 191.475 170.2]);
for i = 1:length(axes_handles)
	hold(axes_handles(i),'on')
end


load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(4).A;
B_spikes = spikes(4).B;
load('/local-data/DA-paper/large-variance-flicker/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
all_spikes = vertcat(all_spikes,spikes(4).A);
B_spikes = vertcat(B_spikes,spikes(4).B);

% A spikes --> firing rate
hash = DataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end

tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 

% plot stimulus
plot(axes_handles(1),tA,mean2(PID),'k')
set(axes_handles(1),'XLim',[10 60],'XTickLabel',{})
ylabel(axes_handles(1),'Stimulus (V)')

% plot response
plot(axes_handles(2),tA,mean2(fA),'k')
set(axes_handles(2),'XLim',[10 60],'XTickLabel',{})
ylabel(axes_handles(2),'Firing Rate (Hz)')

% extract LN model
[K, ~, filtertime_full] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);
K = K/max(K);
fp = convolve(tA,mean2(PID),K,filtertime);
R = mean2(fA);
temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;
hold(axes_handles(2),'on')
clear p
p.A = 57.2707;
p.k = 23.7470;
p.n = 2.9373;
fp_hill = hill(p,fp);

% plot filter + nonlinearity
plot(axes_handles(4),filtertime,K,'r')
plot(axes_handles(5),[0:max(max(fA))],hill(p,[0:max(max(fA))]),'r')
xlabel(axes_handles(4),'Lag (s)')
ylabel(axes_handles(4),'Filter (norm)')
xlabel(axes_handles(5),'Filter Output (Hz)')
ylabel(axes_handles(5),'Nonlinearity (Hz)')

% plot prediction and prediction quality
l=plot(axes_handles(2),tA,fp_hill,'r');
r2 = rsquare(fp_hill,mean2(fA));
legend(l,strcat('r^2=',oval(r2)))


% plot gain
plot(axes_handles(3),tA,mean2(fA)./fp_hill,'r')
ylabel(axes_handles(3),'Gain')
set(axes_handles(3),'XLim',[10 60],'YLim',[0 2])
xlabel(axes_handles(3),'Time (s)')


% gain analysis -- linear model
ph = []; ph(3:4) = axes_handles(6:7);
GainAnalysisWrapper(mean2(fA),fp,mean2(PID),tA,0.4290,ph);

% gain analysis -- LN model
ph = []; ph(3:4) = axes_handles(8:9);
GainAnalysisWrapper(mean2(fA),fp_hill,mean2(PID),tA,0.4290,ph);

% gain aanlysis -- DA model
clear p
p.tau_z = 127.2500;
p.tau_y = 23.8316;
p.  n_y = 2;
p.  n_z = 2;
p.    A = 729.0620;
p.    B = 13.8476;
p.    C = 0.5972;
p.   s0 = -0.1682;
[fp_DA,~,~,Ky,Kz] = DAModelv2(mean2(PID),p);
ph = []; ph(3:4) = axes_handles(10:11);
GainAnalysisWrapper(mean2(fA),fp_DA,mean2(PID),tA,0.4290,ph);

% equalise axes
set(axes_handles(7) ,'XLim',[0.1 10],'YLim',[0 2.5])
set(axes_handles(9) ,'XLim',[0.1 10],'YLim',[0 2.5])
set(axes_handles(11),'XLim',[0.1 10],'YLim',[0 2.5])

% fix some labels
ylabel(axes_handles(7),'Gain')
ylabel(axes_handles(9),'')
ylabel(axes_handles(11),'')
ylabel(axes_handles(6),'Neuron Response (Hz)')
ylabel(axes_handles(8),'')
ylabel(axes_handles(10),'')

xlabel(axes_handles(6),'Linear Prediction (Hz)')
xlabel(axes_handles(8),'LN Prediction (Hz)')
xlabel(axes_handles(10),'DA Prediction (Hz)')

PrettyFig('plw=1.5;','lw=1.5;','fs=14;')


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



