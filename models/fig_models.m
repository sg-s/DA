% fig_models.m
% in which we try to explain observed phenomenon with some models

% created by Srinivas Gorur-Shandilya at 7:11 , 20 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
pHeader;

figure('outerposition',[0 0 1400 900],'PaperUnits','points','PaperSize',[1400 900]); hold on
clear axes_handles
for i = 1:12
	axes_handles(i) = subplot(4,3,i); hold on
end


% ##    ##    ###    ######## ##     ## ########     ###    ##       
% ###   ##   ## ##      ##    ##     ## ##     ##   ## ##   ##       
% ####  ##  ##   ##     ##    ##     ## ##     ##  ##   ##  ##       
% ## ## ## ##     ##    ##    ##     ## ########  ##     ## ##       
% ##  #### #########    ##    ##     ## ##   ##   ######### ##       
% ##   ### ##     ##    ##    ##     ## ##    ##  ##     ## ##       
% ##    ## ##     ##    ##     #######  ##     ## ##     ## ######## 

%  ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
% ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
% ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%  ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%       ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
% ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%  ######     ##    #### ##     ##  #######  ########  #######   ######  


% first we fit the naturalistic stimulus
load('/local-data/DA-paper/natural-flickering/mahmut-raw/2014_07_11_EA_natflick_non_period_CFM_1_ab3_1_1_all.mat')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;

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

% remove the baseline from the PID, and remember the error
PID_baseline = mean(mean(PID(1:5e3,:)));
PID = PID - PID_baseline;


clear data
data.response = mean(fA,2);
data.stimulus = mean(PID,2);
data.response(1:1e3) = NaN;
time = 1e-3*(1:length(data.response));

% first show the data
plot(axes_handles(4),time,data.response,'k')
plot(axes_handles(5),time,data.response,'Color',[0 0 0 0.5])
plot(axes_handles(6),time,data.response,'Color',[0 0 0 0.5])



% show DA model 
clear p
p.   s0 = 3.8594e-04;
p.  n_z = 2;
p.tau_z = 150.4629;
p.  n_y = 2;
p.tau_y = 26.7402;
p.    C = 0.5438;
p.    A = 163.4063;
p.    B = 2.4775;

fp = DAModelv2(data.stimulus,p);
plot(axes_handles(5),time,fp,'r')

% simple receptor model
clear p
p.    r_b = 0.0010;
p.    r_d = 0.0096;
p.theta_b = 2.3984;
p.theta_d = 0.0976;
p. hill_A = 134.5938;
p. hill_K = 0.1377;

fp = simpleReceptorModelv2(data.stimulus,p);
plot(axes_handles(6),time,fp,'r')


% ##     ## ########    ###    ##    ## 
% ###   ### ##         ## ##   ###   ## 
% #### #### ##        ##   ##  ####  ## 
% ## ### ## ######   ##     ## ## ## ## 
% ##     ## ##       ######### ##  #### 
% ##     ## ##       ##     ## ##   ### 
% ##     ## ######## ##     ## ##    ## 

%  ######  ##     ## #### ######## ######## ######## ########  
% ##    ## ##     ##  ##  ##          ##    ##       ##     ## 
% ##       ##     ##  ##  ##          ##    ##       ##     ## 
%  ######  #########  ##  ######      ##    ######   ##     ## 
%       ## ##     ##  ##  ##          ##    ##       ##     ## 
% ##    ## ##     ##  ##  ##          ##    ##       ##     ## 
%  ######  ##     ## #### ##          ##    ######## ########  

%  ######      ###    ##     ##  ######   ######  ####    ###    ##    ##  ######  
% ##    ##    ## ##   ##     ## ##    ## ##    ##  ##    ## ##   ###   ## ##    ## 
% ##         ##   ##  ##     ## ##       ##        ##   ##   ##  ####  ## ##       
% ##   #### ##     ## ##     ##  ######   ######   ##  ##     ## ## ## ##  ######  
% ##    ##  ######### ##     ##       ##       ##  ##  ######### ##  ####       ## 
% ##    ##  ##     ## ##     ## ##    ## ##    ##  ##  ##     ## ##   ### ##    ## 
%  ######   ##     ##  #######   ######   ######  #### ##     ## ##    ##  ######  


[PID, ~, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);


AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)));
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% extract filters and find gain
a = 35e3; z = 55e3;
[K,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

% some core variables
dt = 1e-3;
c = parula(max(paradigm)+1); % colour scheme



% show gain changes for all paradigms -- average over neurons 
ss = 100;
all_x = 0:0.1:2;
axes(axes_handles(7)), hold(axes_handles(7),'on')
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(fA_pred(a:z,paradigm == i),2);
	s = nanmean(PID(a:z,paradigm == i),2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	[~,orn_io_data(i)] = plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end


% fit to all data
clear data
for i = 1:max(paradigm)
	data(i).stimulus = nanmean(PID(a:z,paradigm==i),2);
	data(i).response = nanmean(fA(a:z,paradigm==i),2);
	data(i).response(1:1e3) = NaN;
end



% show DA model i/o curves

clear p
p.   s0 = -0.1373;
p.  n_z = 2;
p.tau_z = 127.8711;
p.  n_y = 2;
p.tau_y = 23.0801;
p.    C = 2.3438e-04;
p.    A = 378.7657;
p.    B = 11.9560;

axes(axes_handles(8))
for i = 1:max(paradigm)
	stim = nanmean(PID(a:z,paradigm==i),2);
	fp = DAModelv2(stim,p);
	fp(1:1e3) = NaN;
	x = nanmean(fA_pred(a:z,paradigm == i),2);
	s = nanmean(PID(a:z,paradigm == i),2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	[~,d] = plotPieceWiseLinear(x,fp,'nbins',50,'make_plot',false);
	errorShade(d.x(1:45),d.y(1:45),d.ye(1:45),'Color',c(i,:));
end

% show receptor model i/o curves
clear p
p.    r_b = 0.0032;
p.    r_d = 0.0211;
p.theta_b = 1.3984;
p.theta_d = 0.0976;
p. hill_A = 136.7032;
p. hill_K = 0.1331;


axes(axes_handles(9))
for i = 1:max(paradigm)
	stim = nanmean(PID(a:z,paradigm==i),2);
	fp = simpleReceptorModelv2(stim,p);
	fp(1:1e3) = NaN;
	x = nanmean(fA_pred(a:z,paradigm == i),2);
	s = nanmean(PID(a:z,paradigm == i),2);
	x = x - nanmean(x);
	x = x + nanmean(nanmean(s));
	[~,d] = plotPieceWiseLinear(x,fp,'nbins',50,'make_plot',false);
	errorShade(d.x(1:45),d.y(1:45),d.ye(1:45),'Color',c(i,:));
end

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
	end


%% Version Info
%
pFooter;


