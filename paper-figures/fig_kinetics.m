% fig_kinetics.m
% 
% created by Srinivas Gorur-Shandilya at 2:14 , 07 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;
dm = dataManager;

figure('outerposition',[0 0 1100 700],'PaperUnits','points','PaperSize',[1100 700]); hold on

% ##    ##    ###    ######## ##     ## ########     ###    ##       
% ###   ##   ## ##      ##    ##     ## ##     ##   ## ##   ##       
% ####  ##  ##   ##     ##    ##     ## ##     ##  ##   ##  ##       
% ## ## ## ##     ##    ##    ##     ## ########  ##     ## ##       
% ##  #### #########    ##    ##     ## ##   ##   ######### ##       
% ##   ### ##     ##    ##    ##     ## ##    ##  ##     ## ##       
% ##    ## ##     ##    ##     #######  ##     ## ##     ## ######## 

%  ######  ######## #### ##     ## ##     ## ##       #### 
% ##    ##    ##     ##  ###   ### ##     ## ##        ##  
% ##          ##     ##  #### #### ##     ## ##        ##  
%  ######     ##     ##  ## ### ## ##     ## ##        ##  
%       ##    ##     ##  ##     ## ##     ## ##        ##  
% ##    ##    ##     ##  ##     ## ##     ## ##        ##  
%  ######     ##    #### ##     ##  #######  ######## #### 


% analyse kinetics of LFP and firing rate during the naturalistic stimulus presentation
load(dm.getPath('aeb361c027b71938021c12a6a12a85cd'),'-mat');

subplot(2,3,4), hold on

c = lines(3);

min_acceptable_corr = .5;
min_acceptable_lag = 2;
for i = 1:length(od)
	S = nanmean(od(i).stimulus,2); 
	R = nanmean(od(i).firing_rate,2);
	X = -nanmean(od(i).LFP,2);

	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,R,1e3,25);
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
	lag(rm_this) = [];
	mean_x(rm_this) = [];

	plotPieceWiseLinear(mean_x,lag,'Color',c(1,:),'nbins',19);


	[lag, mean_x, max_corr] = findLagAndMeanInWindow(S,X,1e3,25);
	rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
	lag(rm_this) = [];
	mean_x(rm_this) = [];

	plotPieceWiseLinear(mean_x,lag,'Color',c(2,:),'nbins',19);
end
xlabel('\mu_{Stimulus} in preceding 1s (V)')
ylabel('Lag (ms)')


% also show an example whiff
a = 2.95e4;
S = S(a:a+500); S = S - min(S); S = S/max(S);
X = X(a:a+500); X = X - min(X); X = X/max(X);
R = R(a:a+500); R = R - min(R); R = R/max(R);
subplot(2,3,1), hold on
title('Naturalistic Stimulus')
plot(S,'k')
plot(X,'Color',c(2,:))
plot(R,'Color',c(1,:))
xlabel('Time since whiff onset (ms)')
legend({'Stimulus','LFP','Firing Rate'})
[~,sp] = max(S);
[~,xp] = max(X);
[~,rp] = max(R);
plot([sp rp],[1.05 1.05],'LineWidth',3,'Color',c(1,:))
plot([sp xp],[1.1 1.1],'LineWidth',3,'Color',c(2,:))


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

clearvars -except being_published dm c
[PID, LFP, fA, paradigm, orn, ~, AllControlParadigms, paradigm_hashes] = consolidateData(dm.getPath('bf79dfd769a97089e42beb0660174e84'),1);

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
LFP(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% band pass all the LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = filtered_LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,filtered_LFP(:,i));
	filtered_LFP(:,i) = filtered_LFP(:,i)*10; % to get the units right, now in mV
end

% some core variables
dt = 1e-3;
lag_LFP = NaN(1,width(PID));
lag_fA = NaN(1,width(PID));
max_corr_LFP = NaN(1,width(PID));
max_corr_fA = NaN(1,width(PID));

subplot(2,3,2), hold on
title('Changing Stimulus Mean')


for i = 1:width(PID)
	s = PID(25e3:45e3,i)-mean(PID(25e3:45e3,i)); s = s/std(s);
	r = fA(25e3:45e3,i)-mean(fA(25e3:45e3,i)); r = r/std(r);
	x = filtered_LFP(25e3:45e3,i)-mean(filtered_LFP(25e3:45e3,i)); x = x/std(x); x = -x;

	[temp,lags] = xcorr(r,s); temp = temp/20e3;
	[max_corr_fA(i),lag_fA(i)] = max(temp);

	if i == 1
		plot(s(1.42e4:1.5e4),'k')
		plot(x(1.42e4:1.5e4),'Color',c(2,:))
		plot(r(1.42e4:1.5e4),'Color',c(1,:))
	end
	% if i == width(PID)
	% 	plot(lags,temp/max(temp),'r--')
	% end

	[temp,lags] = xcorr(x,s);  temp = temp/20e3;
	[max_corr_LFP(i),lag_LFP(i)] = max(temp);

	% if i == 1
	% 	plot(lags,temp/max(temp),'b-')
	% end
	% if i == width(PID)
	% 	plot(lags,temp/max(temp),'b--')
	% end
end
xlabel('Time (ms)')
% ylabel('Cross correlation (norm)')

lag_LFP = lag_LFP - 20e3;
lag_fA = lag_fA - 20e3;
mean_stim =  nanmean(PID(25e3:45e3,:));  
lag_fA(lag_fA<0) = NaN;
lag_LFP(lag_LFP<0) = NaN;

subplot(2,3,5), hold on
plot(mean_stim,lag_LFP,'+','Color',c(2,:))
plot(mean_stim,lag_fA,'+','Color',c(1,:))
xlabel('\mu_{Stimulus} (V)')
ylabel('Lag (ms)')



%  ######   #######  ##    ## ######## ########     ###     ######  ######## 
% ##    ## ##     ## ###   ##    ##    ##     ##   ## ##   ##    ##    ##    
% ##       ##     ## ####  ##    ##    ##     ##  ##   ##  ##          ##    
% ##       ##     ## ## ## ##    ##    ########  ##     ##  ######     ##    
% ##       ##     ## ##  ####    ##    ##   ##   #########       ##    ##    
% ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##    ##    ##    
%  ######   #######  ##    ##    ##    ##     ## ##     ##  ######     ##    

clearvars -except being_published dm c
[PID, LFP, fA, paradigm, orn] = consolidateData(dm.getPath('7955d1ed77512dfe3452b39d71a50e1b'),1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 
% bandpass to remove spikes and slow fluctuations
for i = 1:width(LFP)
	a = find(~isnan(LFP(:,i)),1,'first');
	z = find(~isnan(LFP(:,i)),1,'last');
	LFP(a:z,i) = bandPass(LFP(a:z,i),1000,10)*10; % now in mV
end

% reshape the LFP signals
block_length = 1e4;
reshaped_LFP = LFP(global_start:end-1e4-1,1:width(PID));
reshaped_LFP = reshape(reshaped_LFP,block_length,width(reshaped_LFP)*length(reshaped_LFP)/block_length);

% also reshape the PID
reshaped_PID = PID(global_start:end-1e4-1,1:width(PID));
reshaped_PID = reshape(reshaped_PID,block_length,width(reshaped_PID)*length(reshaped_PID)/block_length);

% reshape the firing rate signals
reshaped_fA = fA(global_start:end-1e4-1,1:width(PID));
reshaped_fA = reshape(reshaped_fA,block_length,width(reshaped_fA)*length(reshaped_fA)/block_length);


% also reshape the orn ID
reshaped_orn = repmat(orn,length(global_start:length(PID)-1e4-1)/block_length,1);
reshaped_orn = reshaped_orn(:);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];
reshaped_orn(rm_this) = [];



lag_LFP = NaN(2,width(reshaped_PID));
lag_fA = NaN(2,width(reshaped_PID));
max_corr_LFP = NaN(2,width(reshaped_PID));
max_corr_fA = NaN(2,width(reshaped_PID));


subplot(2,3,3), hold on
title('Changing Stimulus variance')
for i = 1:width(reshaped_PID)
	s = reshaped_PID(1:5e3,i); s = s - mean(s); s = s/std(s);
	r = reshaped_fA(1:5e3,i); r = r - mean(r); r = r/std(r);
	x = reshaped_LFP(1:5e3,i); x = x - mean(x); x = -x/std(x); 

	[temp,lags] = xcorr(r,s); temp = temp/5e3;
	[max_corr_fA(1,i),lag_fA(1,i)] = max(temp);

	if i == 2
		plot(s(2223:2500),'k')
		plot(x(2223:2500),'Color',c(2,:))
		plot(r(2223:2500),'Color',c(1,:))
	end

	[temp,lags] = xcorr(x,s);  temp = temp/5e3;
	[max_corr_LFP(1,i),lag_LFP(1,i)] = max(temp);

	% if i == 2
	% 	plot(lags,temp/max(temp),'r-')
	% end

	s = reshaped_PID(5e3:end,i); s = s - mean(s); s = s/std(s);
	r = reshaped_fA(5e3:end,i); r = r - mean(r); r = r/std(r);
	x = reshaped_LFP(5e3:end,i); x = x - mean(x); x = -x/std(x); 

	[temp,lags] = xcorr(r,s); temp = temp/5e3;
	[max_corr_fA(2,i),lag_fA(2,i)] = max(temp);

	% if i == 2
	% 	plot(lags,temp/max(temp),'b--')
	% end

	[temp,lags] = xcorr(x,s);  temp = temp/5e3;
	[max_corr_LFP(2,i),lag_LFP(2,i)] = max(temp);

	% if i == 2
	% 	plot(lags,temp/max(temp),'b-')
	% end

end

xlabel('Lag (ms)')

sigma_stim = [std(reshaped_PID(1e3:5e3,:)); std(reshaped_PID(6e3:9e3,:))];

lag_fA = lag_fA - 5e3;
lag_LFP = lag_LFP - 5e3;
lag_LFP(lag_LFP<0) = NaN;
lag_fA(lag_fA<0) = NaN;



subplot(2,3,6), hold on
plot(sigma_stim(1,:),lag_LFP(1,:),'+','Color',c(2,:))
plot(sigma_stim(2,:),lag_LFP(2,:),'+','Color',c(2,:))
plot(sigma_stim(1,:),lag_fA(1,:),'+','Color',c(1,:))
plot(sigma_stim(2,:),lag_fA(2,:),'+','Color',c(1,:))
set(gca,'YLim',[0 200])
xlabel('\sigma_{Stimulus} (V)')
ylabel('Lag (ms)')

legend('boxoff')
prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


