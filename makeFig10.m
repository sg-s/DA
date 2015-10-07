% makeFig10.m
% makes figure 10 of the paper, showing fits of models to the data.
% 
% created by Srinivas Gorur-Shandilya at 9:58 , 07 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


% add homebrew path
path1 = getenv('PATH');
if isempty(strfind(path1,[':/usr/local/bin']))
    path1 = [path1 ':/usr/local/bin'];
end
setenv('PATH', path1);

% this code determines if this function is being called
calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
		unix(['tag -a publish-failed ',which(mfilename)]);
		unix(['tag -r published ',which(mfilename)]);
	end
end
tic

figure('outerposition',[0 0 1400 750],'PaperUnits','points','PaperSize',[1400 750]); hold on
for i = 1:8
	axes_handles(i) = subplot(2,4,i); hold on
end

% get the insets
inset(1) = axes();
set(inset(1),'Position',[.135 .795 .06 .12],'box','on','XTickLabel',{},'YTickLabel',{})

%       ##    ##    ###    ######## ##     ## ########     ###    ##       
%       ###   ##   ## ##      ##    ##     ## ##     ##   ## ##   ##       
%       ####  ##  ##   ##     ##    ##     ## ##     ##  ##   ##  ##       
%       ## ## ## ##     ##    ##    ##     ## ########  ##     ## ##       
%       ##  #### #########    ##    ##     ## ##   ##   ######### ##       
%       ##   ### ##     ##    ##    ##     ## ##    ##  ##     ## ##       
%       ##    ## ##     ##    ##     #######  ##     ## ##     ## ######## 

%        ######  ######## #### ##     ## ##     ## ##       #### 
%       ##    ##    ##     ##  ###   ### ##     ## ##        ##  
%       ##          ##     ##  #### #### ##     ## ##        ##  
%        ######     ##     ##  ## ### ## ##     ## ##        ##  
%             ##    ##     ##  ##     ## ##     ## ##        ##  
%       ##    ##    ##     ##  ##     ## ##     ## ##        ##  
%        ######     ##    #### ##     ##  #######  ######## #### 


load('/local-data/DA-paper/natural-flickering/mahmut-raw/2014_07_11_EA_natflick_non_period_CFM_1_ab3_1_1_all.mat')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;
B_spikes = spikes(2).B;


% A spikes --> firing rate
hash = dataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end

% B spikes --> firing rate
hash = dataHash(full(B_spikes));
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

% fit DA model 
clear p
p.tau_z =  151.0000;
p.tau_y =  26.7063;
p.  n_y =  2;
p.  n_z =  2;
p.    A =  163.2887;
p.    B =  2.4724;
p.    C =  0.5454;
p.   s0 =  -0.1565;
fp_DA = DAModelv2(mean2(PID),p);

% fit a linear model
[K, filtertime_full] = fitFilter2Data(mean2(PID),mean2(fA),'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);
fp_K = convolve(tA,mean2(PID),K,filtertime);


% make the output analysis
shat = computeSmoothedStimulus(mean2(PID),500);
shat = shat-min(shat);
shat = shat/max(shat);
shat = 1+ceil(shat*99);
shat(isnan(shat)) = 1;

ss = 1;
cc = parula(100);
c= cc(shat,:);
R = mean2(fA);
scatter(axes_handles(1),fp_DA(1:ss:end),R(1:ss:end),[],c(1:ss:end,:),'filled')
set(axes_handles(1),'XLim',[0 115],'YLim',[0 115])
xlabel(axes_handles(1),'DA Model Prediction (Hz)')
ylabel(axes_handles(1),'ORN Response (Hz)')

scatter(inset(1),fp_K(1:ss:end),R(1:ss:end),[],c(1:ss:end,:),'filled')
set(inset(1),'box','on','XLim',[min(fp_K) max(fp_K)],'YLim',[min(R) max(R)],'XTickLabel',{},'YTickLabel',{})

% quantify the gain for the linear prediction
% first account for some trivial scaling
temp =fit(fp_K(~(isnan(fp_K) | isnan(R))),R(~(isnan(fp_K) | isnan(R))),'poly1');
fp_K = fp_K*temp.p1;
fp_K = fp_K+temp.p2;

shat = computeSmoothedStimulus(mean2(PID),500);
[whiff_starts,whiff_ends] = computeOnsOffs(R>10);
mean_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err = gain;
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	ff=fit(fp_K(whiff_starts(i):whiff_ends(i)),R(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
rm_this = (abs(gain_err./gain)) > .5; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
mean_stim(rm_this) = [];
l(1) = plot(axes_handles(5),mean_stim,gain,'k+');
xlabel(axes_handles(5),'Mean Stimulus in last 500ms')
ylabel(axes_handles(5),'Relative Gain')

% now do the same for the DA model
gain = NaN*whiff_ends;
gain_err = gain;
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	ff=fit(fp_DA(whiff_starts(i):whiff_ends(i)),R(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
rm_this = (abs(gain_err./gain)) > .5; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
mean_stim(rm_this) = [];
l(2) = plot(axes_handles(5),mean_stim,gain,'r+');
legend(l,{'Linear Prediction','DA Model'})


prettyFig('fs=18;','lw=1.5;')

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
disp(dataHash(strcat(mfilename,'.m'),Opt))

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

% tag the file as being published 

if being_published
	unix(['tag -a published ',which(mfilename)]);
	unix(['tag -r publish-failed ',which(mfilename)]);
end
