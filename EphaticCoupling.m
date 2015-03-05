% EphaticCoupling.m
% 
% created by Srinivas Gorur-Shandilya at 10:15 , 05 March 2015. Contact me at http://srinivas.gs/contact/
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


%% Ephatic Coupling
% In this document, we investigate the ephatic coupling between the A and B neuron as shown in Chih-Ying's  
% <http://www.nature.com/nature/journal/v492/n7427/full/nature11712.html paper>

%% Data
% We use the very clean data from the the Large Variance Flicker experiment: 

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


%        ########  ########  ######  ########   #######  ##    ##  ######  ######## 
%        ##     ## ##       ##    ## ##     ## ##     ## ###   ## ##    ## ##       
%        ##     ## ##       ##       ##     ## ##     ## ####  ## ##       ##       
%        ########  ######    ######  ########  ##     ## ## ## ##  ######  ######   
%        ##   ##   ##             ## ##        ##     ## ##  ####       ## ##       
%        ##    ##  ##       ##    ## ##        ##     ## ##   ### ##    ## ##       
%        ##     ## ########  ######  ##         #######  ##    ##  ######  ######## 




figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
subplot(2,9,10:14), hold on
plot(tA,mean2(fA),'k')
set(gca,'XLim',[10 60])
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


%          ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
%         ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
%         ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%          ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%               ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
%         ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%          ######     ##    #### ##     ##  #######  ########  #######   ######  


subplot(2,9,1:5), hold on
plot(tA,mean2(PID),'k')
set(gca,'XLim',[10 60])
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
set(gca,'XLim',[10 60])
xlabel('Time (s)')
ylabel('Firing Rate (A) (Hz)')
set(gca,'YLim',[0 70])

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
set(gca,'XLim',[10 60])
xlabel('Time (s)')
ylabel('Firing Rate (B) (Hz)')
set(gca,'YLim',[0 70])

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
% Let's look a raster plot of the spikes:

figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
raster2(all_spikes,B_spikes)

xlabel('Time (s)')
ylabel('Trial #')
PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% It looks like there is something weird about the 2nd and 3rd trials. We will skip that from the analysis. 
skip  = [2 3];

%%
% First, we simply measure the coefficient of determination between the A and B neurons:

disp(rsquare(mean2(fA),mean2(fB)))

%%
% Now, we plot the cross correlation between the A and B neurons: 

dt = mean(diff(tA));
x = NaN(width(fA),2001);
tc = zeros(2001,1);
for i = setdiff(1:width(fA),skip)
	temp=xcorr(full(fB(:,i))-mean(fB(:,i)),full(fA(:,i))-mean(fA(:,i)));
	temp = temp/max(temp(length(tA)-1e3:length(tA)+1e3));
	tt = dt*(1:length(temp));
	tt = tt -mean(tt);
	tc = tt(length(tA)-1e3:length(tA)+1e3);
	x(i,:)= temp(length(tA)-1e3:length(tA)+1e3);
end


figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
plot(tc,x)
xlabel('Lag (s) A \rightarrow B')
ylabel('Cross Correlation (norm)')

PrettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% This is something weird. Other than one trial (which is the first trial of the 2nd neuron), every trial shows a positive cross correlation around zero. Based on ephatic coupling, would we not expect that the value around zero to be most negative? 

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