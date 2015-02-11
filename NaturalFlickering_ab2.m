% Mahmut_Data_Analyis.m
% this is a complete rewrite on  4:06 , 10 February 2015. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
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

load('/local-data/DA-paper/natural-flickering/mahmut-raw/2014_07_11_EA_natflick_non_period_CFM_1_ab2_1_1_all.mat')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;
B_spikes = spikes(2).B;


% A spikes --> firing rate
hash = DataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end

% B spikes --> firing rate
hash = DataHash(full(B_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fB = spiketimes2f(B_spikes,time);
	cache(hash,fB);
else
	fB = cached_data;
end

tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 

% trash first trial
fA(:,1) = [];
PID(:,1) = [];



figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
subplot(2,9,10:14), hold on
plot(tA,mean2(fA),'k')
set(gca,'XLim',[0 70])
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(2,9,15:16), hold on
[r2,s] = rsquare(fA);
imagescnan(r2)
caxis([0 1])
axis image
colorbar
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))

subplot(2,9,17:18), hold on
imagescnan(s)
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))


subplot(2,9,1:5), hold on
plot(tA,mean2(PID),'k')
set(gca,'YScale','log')
set(gca,'XLim',[0 70])
xlabel('Time (s)')
ylabel('Odor Concentration (V)')

subplot(2,9,6:7), hold on
[r2,s] = rsquare(PID);
imagescnan(r2)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))

subplot(2,9,8:9), hold on
imagescnan(s)
colorbar
axis image
axis off
title(strcat('mean slope = ',oval(mean(s(~isnan(s))),2)))

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%%
% For completeness, here is a comparison of the A neuron's response to that of the B neuron's. 


figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
subplot(2,8,1:6), hold on
plot(tA,mean2(fA),'k')
set(gca,'XLim',[0 70])
xlabel('Time (s)')
ylabel('Firing Rate (A) (Hz)')
set(gca,'YLim',[0 200])

subplot(2,8,7:8), hold on
hash = DataHash(fA);
cached_data = cache(hash);
if isempty(cached_data)
	r2 = rsquare(fA);
	cache(hash,r2);
else
	r2 = cached_data;
end
imagescnan(r2)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))


subplot(2,8,9:14), hold on
plot(tA,mean2(fB),'k')
set(gca,'XLim',[0 70])
xlabel('Time (s)')
ylabel('Firing Rate (B) (Hz)')
set(gca,'YLim',[0 200])

subplot(2,8,15:16), hold on
hash = DataHash(fB);
cached_data = cache(hash);
if isempty(cached_data)
	r2 = rsquare(fB);
	cache(hash,r2);
else
	r2 = cached_data;
end
imagescnan(r2)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('mean r^2 = ',oval(mean(r2(~isnan(r2))),2)))


PrettyFig();

if being_published
	snapnow
	delete(gcf)
end



%%
% In this document, we generate odor stimuli flickers over a large range, like the "natural" stimuli, but never goes to zero, so that the neuron should never silence (allowing us to accurately follow its response, and reasonably estimate instantaneous gain). 

%%
% The first figure shows the odor stimulus (ethyl acetate) presented to two ab3A neurons. Each neuron was recorded from ten times. The colormaps on the right show the coefficient of determination between pairwise trials. 


%%
% The following figure shows the stimulus distribution and the autocorrelation functions of the stimulus and the response. 


figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1), hold on
for i = 1:width(PID)
	[y,x] = histcounts(PID(:,i),300);x(1) = [];
	plot(x,y)
end

set(gca,'YScale','log','XScale','log')
xlabel('PID (V)')
ylabel('Count')
title('Stimulus Histogram')

subplot(1,3,2), hold on
a = [];
for i = 1:width(PID)
	[~,~,c,temp]=FindCorrelationTime(PID(:,i));
	a = [a temp];
	l=plot(1:3e3,c(1:3e3));
end

xlabel('Lag (ms)')
ylabel('Autocorrelation')
legend(l,strcat('\tau=',oval(a),'ms'))
set(gca,'XScale','log')
title('Stimulus')

subplot(1,3,3), hold on
a = [];
for i = 1:width(fA)
	[~,~,c,temp]=FindCorrelationTime(fA(:,i));
	a = [a temp];
	l=plot(1:3e3,c(1:3e3));
end

title('Response')
xlabel('Lag (ms)')
ylabel('Autocorrelation')
legend(l,strcat('\tau=',oval(a),'ms'))
set(gca,'XScale','log')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end



%     ##       #### ##    ## ########    ###    ########         ######## #### ######## 
%     ##        ##  ###   ## ##         ## ##   ##     ##        ##        ##     ##    
%     ##        ##  ####  ## ##        ##   ##  ##     ##        ##        ##     ##    
%     ##        ##  ## ## ## ######   ##     ## ########         ######    ##     ##    
%     ##        ##  ##  #### ##       ######### ##   ##          ##        ##     ##    
%     ##        ##  ##   ### ##       ##     ## ##    ##         ##        ##     ##    
%     ######## #### ##    ## ######## ##     ## ##     ##        ##       ####    ##    


%% Linear Fit
% In this section we fit a linear kernel, first to the stimulus, and then to the log of the stimulus using standard methods. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[900 500]); hold on
subplot(1,2,1), hold on
plot([-.1 1],[0 0 ],'k--')
[K, ~, filtertime_full] = FindBestFilter(mean2(PID),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

plot(filtertime,K,'r')
ylabel('Filter Amplitude')
xlabel('Filter Lag (s)')
title('Filter from stimulus')

subplot(1,2,2), hold on
plot([-.1 1],[0 0 ],'k--')
[Klog, ~, filtertime_full] = FindBestFilter(log(mean2(PID)),mean2(fA),[],'regmax=1;','regmin=1;','filter_length=1999;','offset=500;');
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
Klog = interp1(filtertime_full,Klog,filtertime);

plot(filtertime,Klog,'r')
ylabel('Filter Amplitude')
xlabel('Filter Lag (s)')
title('Filter from log stimulus')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%%
% We will use this filter to make a prediction of the response:

% convolve with filter to make prediction
fp = convolve(tA,mean2(PID),K,filtertime);
fp_log = convolve(tA,log(mean2(PID)),Klog,filtertime);

% correct for trivial scaling
R = mean2(fA);
temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;
temp =fit(fp_log(~(isnan(fp_log) | isnan(R))),R(~(isnan(fp_log) | isnan(R))),'poly1');
fp_log = fp_log*temp.p1;
fp_log = fp_log+temp.p2;


figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on
subplot(2,4,1:3), hold on
plot(tA,mean2(fA),'k')
l=plot(tA,fp,'r');
r2 = rsquare(fp,mean2(fA));
legend(l,strcat('r^2=',oval(r2,2)))
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(2,4,4), hold on
plot(filtertime,K,'r')
xlabel('Filter Lag (s)')
ylabel('Filter (norm.)')
title('Filter from stimulus')

subplot(2,4,5:7), hold on
plot(tA,mean2(fA),'k')
l=plot(tA,fp_log,'r');
r2 = rsquare(fp_log,mean2(fA));
legend(l,strcat('r^2=',oval(r2,2)))
xlabel('Time (s)')
ylabel('Firing Rate (Hz)')

subplot(2,4,8), hold on
plot(filtertime,Klog,'r')
xlabel('Filter Lag (s)')
ylabel('Filter (norm.)')
title('Filter from log stimulus')

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end

%             #######  ##     ## ######## ########  ##     ## ######## 
%            ##     ## ##     ##    ##    ##     ## ##     ##    ##    
%            ##     ## ##     ##    ##    ##     ## ##     ##    ##    
%            ##     ## ##     ##    ##    ########  ##     ##    ##    
%            ##     ## ##     ##    ##    ##        ##     ##    ##    
%            ##     ## ##     ##    ##    ##        ##     ##    ##    
%             #######   #######     ##    ##         #######     ##    
           
%               ###    ##    ##    ###    ##       ##    ##  ######  ####  ######  
%              ## ##   ###   ##   ## ##   ##        ##  ##  ##    ##  ##  ##    ## 
%             ##   ##  ####  ##  ##   ##  ##         ####   ##        ##  ##       
%            ##     ## ## ## ## ##     ## ##          ##     ######   ##   ######  
%            ######### ##  #### ######### ##          ##          ##  ##        ## 
%            ##     ## ##   ### ##     ## ##          ##    ##    ##  ##  ##    ## 
%            ##     ## ##    ## ##     ## ########    ##     ######  ####  ######  


%% Output Analysis
% In this section, we analyse the prediction of the linear kernel in some more detail. The following figure shows a plot of the linear predictions vs. the actual response. 

ss = 10;
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(fp(1:ss:end),R(1:ss:end),'.','MarkerSize',12,'MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor',[0.9 0.9 0.9])
xlabel('Linear Prediction (Hz)')
ylabel('Actual response (Hz)')
clear p l
p.A = 149.8905;
p.k = 43.2554;
p.n = 2.1188;
fp_LN = hill(p,fp);
l = plot(1:300,hill(p,1:300),'r');
r2 = strcat('r^2=',oval(rsquare(fp_LN,R)));
legend(l,r2,'Location','southeast');

subplot(1,2,2), hold on
plot(fp_log(1:ss:end),R(1:ss:end),'.','MarkerSize',12,'MarkerFaceColor',[0.9 0.9 0.9],'MarkerEdgeColor',[0.9 0.9 0.9])
xlabel('Linear Prediction (log stim.) (Hz)')
clear p l
p.A = 218.4697;
p.k = 93.7756;
p.n = 1.4129;
fp_log_LN = hill(p,fp_log);
l = plot(1:300,hill(p,1:300),'r');
r2 = strcat('r^2=',oval(rsquare(fp_log_LN,R)));
legend(l,r2,'Location','southeast');

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


%% 
% We clearly see that there are large loops in this space, that cannot be fit by any static nonlinear function. These loops are characteristic signatures of a dynamical process that escapes a full description by the linear kernel. 

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


	