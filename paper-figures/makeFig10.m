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

figure('outerposition',[0 0 1600 750],'PaperUnits','points','PaperSize',[1600 750]); hold on
for i = 1:10
	axes_handles(i) = subplot(2,5,i); hold on
end

% get the insets
inset(1) = axes(); % to show natural flickering vs. linear prediction
set(inset(1),'Position',[.2 .6 .05 .12],'box','on','XTickLabel',{},'YTickLabel',{})

inset(2) = axes();
set(inset(2),'Position',[.688 .795 .05 .12],'box','on','XTickLabel',{},'YTickLabel',{})

inset(3) = axes();
set(inset(3),'Position',[.688 .322 .05 .12],'box','on','XTickLabel',{},'YTickLabel',{})

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

% fit a linear model to the orn response
[K, filtertime_full] = fitFilter2Data(mean2(PID),mean2(fA),'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);
fp_K = convolve(tA,mean2(PID),K,filtertime);

% fit another linear model to the DA model
K2 = fitFilter2Data(mean2(PID),fp_DA,'reg',1,'filter_length',1999,'offset',500);
K2 = interp1(filtertime_full,K2,filtertime);
fp_K2 = convolve(tA,mean2(PID),K2,filtertime);

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

scatter(inset(1),fp_K(1:ss:end),R(1:ss:end),[],c(1:ss:end,:),'filled')

scatter(axes_handles(1),fp_K2(1:ss:end),fp_DA(1:ss:end),[],c(1:ss:end,:),'filled')
set(axes_handles(1),'XLim',[min(fp_K2) max(fp_K2)],'YLim',[0 115])
xlabel(axes_handles(1),'Projected Stimulus')
ylabel(axes_handles(1),'DA Model Prediction (Hz)')

set(inset(1),'box','on','XLim',[min(fp_K) max(fp_K)],'YLim',[min(R) max(R)],'XTickLabel',{},'YTickLabel',{})


% quantify the gain for the linear prediction
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
l(1) = plot(axes_handles(6),mean_stim,gain,'k+');
xlabel(axes_handles(6),'Mean Stimulus in last 500ms (V)')
ylabel(axes_handles(6),'Relative Gain')

% now do the same for the DA model
gain = NaN*whiff_ends;
gain_err = gain;
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	ff=fit(fp_K2(whiff_starts(i):whiff_ends(i)),fp_DA(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
rm_this = (abs(gain_err./gain)) > .5; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
mean_stim(rm_this) = [];
l(2) = plot(axes_handles(6),mean_stim,gain,'r+');
xlabel(axes_handles(6),'Mean Stimulus in last 500ms (V')
ylabel(axes_handles(6),'Gain (Hz/V)')
set(axes_handles(6),'YLim',[0 165])
legend(l,{'ORN Response','DA Model'})

%       ##      ## ######## ########  ######## ########          
%       ##  ##  ## ##       ##     ## ##       ##     ##         
%       ##  ##  ## ##       ##     ## ##       ##     ##         
%       ##  ##  ## ######   ########  ######   ########  ####### 
%       ##  ##  ## ##       ##     ## ##       ##   ##           
%       ##  ##  ## ##       ##     ## ##       ##    ##          
%        ###  ###  ######## ########  ######## ##     ##         
      
%       ######## ########  ######  ##     ## ##    ## ######## ########  
%       ##       ##       ##    ## ##     ## ###   ## ##       ##     ## 
%       ##       ##       ##       ##     ## ####  ## ##       ##     ## 
%       ######   ######   ##       ######### ## ## ## ######   ########  
%       ##       ##       ##       ##     ## ##  #### ##       ##   ##   
%       ##       ##       ##    ## ##     ## ##   ### ##       ##    ##  
%       ##       ########  ######  ##     ## ##    ## ######## ##     ## 

clearvars -except axes_handles inset being_published 
load('../data/MeanShiftedGaussians.mat')
load('../data/MSG_per_neuron.mat','MSG_data')
load('.cache/hill_fits_MSG.mat','p')
c = parula(9);

% show fits of Hill functions to each cloud
all_x = [0 5];
all_x = linspace(all_x(1),all_x(end),100);
for i = 1:8
	s = ([MSG_data(i,:).stim]);
	if ~isvector(s)
		s = mean2(s);
	end
	filtertime = (-200:800)*1e-3;
	K = reshape([MSG_data(i,:).K],1001,length([MSG_data(i,:).K])/1001);
	if ~isvector(K)
		K = mean2(K);
	end 
	x = convolve(MSG_data(1,1).time,s,K,filtertime) + mean(s);
	plot(axes_handles(2),all_x(all_x<max(x)),hill(all_x(all_x<max(x)),p(i)),'Color',c(i,:))
end

xlabel(axes_handles(2),'Projected Stimulus (V)')
ylabel(axes_handles(2),'Output of nonlinearity (Hz)')

% plot the gain of the neuron vs mean stimulus
% compute gain changes on a per-neuron basis
gain = NaN(8,13);
gain_LN = NaN(8,13);
mean_stim = NaN(8,13);
gain_err = NaN(8,13);
gain_err_LN = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			y = MSG_data(i,j).resp; % average over all neurons 
			s = ([MSG_data(i,:).stim]);
			if ~isvector(s)
				s = mean2(s);
			end
			if ~isvector(y)
				y = mean2(y);
			end 

			K = reshape([MSG_data(i,j).K],1001,length([MSG_data(i,j).K])/1001);
			x = convolve(MSG_data(1,1).time,s,K,filtertime) + mean(s);
			
			% trim NaNs again
			rm_this = isnan(x) | isnan(y);
			x(rm_this) = [];
			y(rm_this) = [];

			temp=fit(x(:),y(:),'poly1');
			gain(i,j) = temp.p1;
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));

			% make a LN prediction -- average over all hill functions
			K = reshape([MSG_data(i,j).K],1001,length([MSG_data(i,j).K])/1001);
			x = convolve(MSG_data(1,1).time,s,K,filtertime) + mean(s);
			clear y
			for k = 1:8
				y(:,k) = hill(x,p(k));
			end
			y = mean2(y);

			% trim NaNs again
			rm_this = isnan(x) | isnan(y);
			x(rm_this) = [];
			y(rm_this) = [];

			temp=fit(x(:),y(:),'poly1');
			gain_LN(i,j) = temp.p1;

		end
	end	
end

% show gain changes -- gain vs. mean stimulus
for i = 1:8
	plot(axes_handles(7),(mean_stim(i,:)),(gain(i,:)),'+','Color',c(i,:));
end
set(axes_handles(7),'YScale','log','XScale','log')


% fit a DA model to all the data
% clear d
% for j = 2
% 	clear d c
% 	c = 1;
% 	for i = 1:8
% 		if ~isempty(MSG_data(i,j).stim)
% 			stim = ([MSG_data(i,j).stim]);
% 			resp = ([MSG_data(i,j).resp]);
% 			if width(stim) > 1
% 				d(c).stimulus = mean2(stim);
% 				d(c).response = mean2(resp);
% 			else
% 				d(c).stimulus = (stim);
% 				d(c).response = (resp);
% 			end
% 			d(c).response(1:1e3) = NaN;
% 			c = c+1;
			
% 		end
% 	end
% 	p(j) = fitModel2Data(@DAModelv2,d,'nsteps',300,'p0',p(j));
% end
load('../data/DA_Fit_to_MSG.mat')

% make DA prediction everywhere
for i = 1:8
	for j = 1:13
		if ~isempty(MSG_data(i,j).stim)
			stim = ([MSG_data(i,j).stim]);
			if width(stim)>1
				stim = mean2(stim);
			end
			MSG_data(i,j).fp_DA = DAModelv2(stim,p(j));
		end
	end
end

% now compute gains for the DA model fit to this data
gain_DA = NaN(8,13);
for i = 1:8 % iterate over all paradigms 
	for j = 1:13
		if width(MSG_data(i,j).stim) > 1
			s = ([MSG_data(i,:).stim]);
			if ~isvector(s)
				s = mean2(s);
			end
			y = MSG_data(i,j).fp_DA;
			if ~isvector(y)
				y = mean2(y);
			end 

			K = reshape([MSG_data(i,j).K],1001,length([MSG_data(i,j).K])/1001);
			x = convolve(MSG_data(1,1).time,s,K,filtertime) + mean(s);

			% trim NaNs again
			rm_this = isnan(x) | isnan(y);
			x(rm_this) = [];
			y(rm_this) = [];

			temp=fit(x(:),y(:),'poly1');
			gain_DA(i,j) = temp.p1;
			mean_stim(i,j) = mean(mean([MSG_data(i,j).stim]));
		end
	end	
end

plot(axes_handles(7),mean_stim(:),gain_DA(:),'r+')
xlabel(axes_handles(7),'Mean Stimulus (V)')
ylabel(axes_handles(7),'ORN Gain (Hz/V)')

%       ######## #### ##       ######## ######## ########  
%       ##        ##  ##          ##    ##       ##     ## 
%       ##        ##  ##          ##    ##       ##     ## 
%       ######    ##  ##          ##    ######   ########  
%       ##        ##  ##          ##    ##       ##   ##   
%       ##        ##  ##          ##    ##       ##    ##  
%       ##       #### ########    ##    ######## ##     ## 
      
%        ######  ########  ######## ######## ########  ##     ## ########   ######  
%       ##    ## ##     ## ##       ##       ##     ## ##     ## ##     ## ##    ## 
%       ##       ##     ## ##       ##       ##     ## ##     ## ##     ## ##       
%        ######  ########  ######   ######   ##     ## ##     ## ########   ######  
%             ## ##        ##       ##       ##     ## ##     ## ##              ## 
%       ##    ## ##        ##       ##       ##     ## ##     ## ##        ##    ## 
%        ######  ##        ######## ######## ########   #######  ##         ######  



% compute peak locations of all these filters
clear l 
l = zeros(8,1);
peak_loc_K = NaN(8,13);
mean_stim_K = NaN(8,13);
for i = 1:8
	for j = 1:13
		if ~isempty(MSG_data(i,j).K)
			K2 = pFilter(MSG_data(i,j).K(200:end),MSG_data(i,j).p);
			filtertime = 1e-3*(1:length(K2));
			l(i)=plot(axes_handles(3),filtertime,K2,'Color',c(i,:));
			[~,loc] = max(K2);
			peak_loc_K(i,j) = filtertime(loc);
			mean_stim_K(i,j) = mean(mean(MSG_data(i,j).stim));
		end
	end
end

set(axes_handles(3),'XLim',[-.01 .5])
xlabel(axes_handles(3),'Lag (s)')
ylabel(axes_handles(3),'Filter (norm)')
L = paradigm_names;
for i = 1:length(L)
	L{i} = L{i}(strfind(L{i},'-')+1:end);
end
legend(l,L);


% extract filters for each case 
for i = 1:8
	for j = 1:13
		if ~isempty(MSG_data(i,j).stim)
			stim = ([MSG_data(i,j).stim]);
			if width(stim)>1
				stim = mean2(stim);
			end
			MSG_data(i,j).K_DA = fitFilter2Data(stim,MSG_data(i,j).fp_DA,'reg',1);
		end
	end
end

% find filter peaks for each neuron, for each paradigm 
peak_loc_K_DA = NaN(8,13);
for i = 1:8
	for j = 1:13
		if ~isempty(MSG_data(i,j).K_DA)
			K2 = MSG_data(i,j).K_DA;
			filtertime = 1e-3*(1:length(K2));
			[~,peak_loc_K_DA(i,j)] = max(K2);
		end
	end
end

clear l
l(1) = plot(axes_handles(8),nonnans(mean_stim_K(:)),nonnans(peak_loc_K(:)*1e3),'k+');
l(2) = plot(axes_handles(8),mean_stim_K(~isnan(peak_loc_K_DA)),nonnans(peak_loc_K_DA(:)),'r+');

set(axes_handles(8),'YLim',[30 120])
xlabel(axes_handles(8),'Mean Stimulus (V)')
ylabel(axes_handles(8),'Peak Time (ms)')


%       ########    ###     ######  ########     ######      ###    #### ##    ## 
%       ##         ## ##   ##    ##    ##       ##    ##    ## ##    ##  ###   ## 
%       ##        ##   ##  ##          ##       ##         ##   ##   ##  ####  ## 
%       ######   ##     ##  ######     ##       ##   #### ##     ##  ##  ## ## ## 
%       ##       #########       ##    ##       ##    ##  #########  ##  ##  #### 
%       ##       ##     ## ##    ##    ##       ##    ##  ##     ##  ##  ##   ### 
%       ##       ##     ##  ######     ##        ######   ##     ## #### ##    ## 
      
%        ######   #######  ##    ## ######## ########   #######  ##       
%       ##    ## ##     ## ###   ##    ##    ##     ## ##     ## ##       
%       ##       ##     ## ####  ##    ##    ##     ## ##     ## ##       
%       ##       ##     ## ## ## ##    ##    ########  ##     ## ##       
%       ##       ##     ## ##  ####    ##    ##   ##   ##     ## ##       
%       ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##       
%        ######   #######  ##    ##    ##    ##     ##  #######  ######## 

clearvars -except inset axes_handles being_published

% the following code was copied from an old commit-- go over this carefully 

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
hash = dataHash(full(all_spikes));
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

% extract LN model
[K, filtertime_full] = fitFilter2Data(mean2(PID),mean2(fA),'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

fp = convolve(tA,mean2(PID),K,filtertime);
R = mean2(fA);
temp =fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;

clear p
p.A = 57.1534;
p.k = 23.6690;
p.n = 2.9341;
fp_hill = hill(p,fp);


% gain analysis -- LN model
clear ph
ph(4) = axes_handles(4);
ph(3) = inset(2);

hl_min = .1;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);

resp = mean2(fA(10e3:55e3,[3:10 13:20]));
pred = fp_hill(10e3:55e3);
time = 1e-3*(1:length(resp));
stim = mean2(PID(10e3:55e3,[3:10 13:20]));

[p,~,~,~,~,history_lengths]=gainAnalysisWrapper('response',resp,'prediction',pred,'stimulus',stim,'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',.5,'use_cache',true,'engine',@gainAnalysis);

ylabel(axes_handles(4),'Gain cf. LN model')
title(axes_handles(4),'')
title(inset(2),'')
set(inset(2),'XScale','linear','XTickLabel','','YTickLabel','')

set(axes_handles(4),'YLim',[0.5 2])

% thin some of the scatter points
lines = get(inset(2),'Children');
delete(lines(5))

warning off
lines(3).YData = lines(3).YData(1:5:end);
lines(3).XData = lines(3).XData(1:5:end);
lines(4).YData = lines(4).YData(1:5:end);
lines(4).XData = lines(4).XData(1:5:end);
warning on

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
fp_DA = DAModelv2(mean2(PID(10e3:55e3,[3:10 13:20])),p);

% gain analysis -- DA model
clear ph
ph(4) = axes_handles(9);
ph(3) = inset(3);

resp = mean2(fA(10e3:55e3,[3:10 13:20]));
pred = fp_DA;
time = 1e-3*(1:length(resp));
stim = mean2(PID(10e3:55e3,[3:10 13:20]));

[p,~,~,~,~,history_lengths] = gainAnalysisWrapper('response',resp,'prediction',pred,'stimulus',stim,'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',.5,'use_cache',true,'engine',@gainAnalysis);

ylabel(axes_handles(9),'Gain cf. DA model')
title(axes_handles(9),'')
title(inset(3),'')
set(inset(3),'XScale','linear','XTickLabel','','YTickLabel','')

set(axes_handles(9),'YLim',[0.5 2])

% thin some of the scatter points
lines = get(inset(3),'Children');
delete(lines(5))

warning off
lines(3).YData = lines(3).YData(1:5:end);
lines(3).XData = lines(3).XData(1:5:end);
lines(4).YData = lines(4).YData(1:5:end);
lines(4).XData = lines(4).XData(1:5:end);
warning on


prettyFig('fs=18;','lw=1.5;')

inset_font_size = 13;
xlabel(inset(2),'LN Prediction','FontSize',inset_font_size)
ylabel(inset(2),'ORN Response','FontSize',inset_font_size)
xlabel(inset(3),'DA Prediction','FontSize',inset_font_size)
ylabel(inset(3),'ORN Response','FontSize',inset_font_size)

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
