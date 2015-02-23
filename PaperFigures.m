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

%% Figure 2

%% Figure 3

clearvars -except being_published
fig_handle=figure('Units','pixels','outerposition',[66 5 661 871],'Color','w','PaperUnits','points','PaperSize',[661 871],'Color','w','Toolbar','none');
clf(fig_handle);
axes_handles(1)=axes('Units','pixels','Position',[63.825 744.625 574.425 85.1]);
axes_handles(2)=axes('Units','pixels','Position',[63.825 595.7 574.425 127.65]);
axes_handles(3)=axes('Units','pixels','Position',[63.825 489.325 574.425 85.1]);
axes_handles(4)=axes('Units','pixels','Position',[63.825 276.575 191.475 170.2]);
axes_handles(5)=axes('Units','pixels','Position',[382.95 276.575 255.3 170.2]);
axes_handles(6)=axes('Units','pixels','Position',[63.825 42.55 191.475 191.475]);
axes_handles(7)=axes('Units','pixels','Position',[382.95 42.55 212.75 191.475]);
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

% extract Linear model
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

% % plot filter + nonlinearity
% plot(axes_handles(4),filtertime,K,'r')
% plot(axes_handles(5),[0:max(max(fA))],hill(p,[0:max(max(fA))]),'r')
% xlabel(axes_handles(4),'Lag (s)')
% ylabel(axes_handles(4),'Filter (norm)')
% xlabel(axes_handles(5),'Filter Output (Hz)')
% ylabel(axes_handles(5),'Nonlinearity (Hz)')

% plot prediction and prediction quality
l=plot(axes_handles(2),tA,fp,'r');
r2 = rsquare(fp,mean2(fA));
legend(l,strcat('r^2=',oval(r2)))


% plot gain
plot(axes_handles(3),tA,mean2(fA)./fp,'r')
ylabel(axes_handles(3),'Gain')
set(axes_handles(3),'XLim',[10 60],'YLim',[0 2])
xlabel(axes_handles(3),'Time (s)')


% gain analysis -- linear model
ph = []; ph(3:4) = axes_handles(4:5);
GainAnalysisWrapper(mean2(fA),fp,mean2(PID),tA,0.4290,ph);

% plot gain vs preceding stimulus
[x,y] = MakeFig6G(mean2(PID),mean2(fA),fp,400);
rm_this = (isnan(x) | isnan(y));
x(rm_this) = [];
y(rm_this) = [];
ss = 50;
plot(axes_handles(6),x(1:ss:end),y(1:ss:end),'k.')
xlabel(axes_handles(6),'Stimulus in preceding 400ms')
ylabel(axes_handles(6),'Instantaneous gain')


% fix some labels
ylabel(axes_handles(5),'Gain')
xlabel(axes_handles(4),'Linear Prediction (Hz)')
ylabel(axes_handles(4),'Neuron Response (Hz)')

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



