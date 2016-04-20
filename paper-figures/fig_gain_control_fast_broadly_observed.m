% fig_gain_control_fast_broadly_observed
% makes a figure for the paper showing that gain control is fast and is broadly observed


pHeader;

figure('outerposition',[0 0 1600 700],'PaperUnits','points','PaperSize',[1600 700]); hold on

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


load('/local-data/DA-paper/natural-flickering/without-lfp/2014_07_11_EA_natflick_non_period_CFM_1_ab3_1_1_all.mat')
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

% make a linear filter
R = mean(fA,2);
[K, filtertime_full] = fitFilter2Data(mean(PID,2),R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% convolve with filter to make prediction
fp = convolve(tA,mean(PID,2),K,filtertime);

shat = computeSmoothedStimulus(mean(PID,2),500);

% find all excursions (defined as firing rate crossing 10Hz)
[whiff_starts,whiff_ends] = computeOnsOffs(R>10);
mean_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err =  NaN*whiff_ends;
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	ff=fit(fp(whiff_starts(i):whiff_ends(i)),R(whiff_starts(i):whiff_ends(i)),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
rm_this = (abs(gain_err./gain)) > .5; % throw out points where the estimate of gain has a more than 50% error
gain(rm_this) = [];
gain_err(rm_this) = [];
whiff_starts(rm_this) = [];
mean_stim(rm_this) = [];

subplot(2,5,1), hold on
plot(mean_stim,gain,'k+');
set(gca,'XScale','log','YScale','log')
xlabel('\mu_{Stimulus} in preceding 500ms (V)')
ylabel('ORN Gain (Hz)')

% compute rho for various history lengths
tau_h = round(logspace(1,4,50));
rho = computeRhoForGainTimescales(gain,whiff_starts,mean(PID,2),tau_h);
subplot(2,5,6), hold on 
plot(tau_h,rho,'k+')
set(gca,'XScale','log','XTick',[10 1e2 1e3 1e4])
xlabel('History Length (ms)')
ylabel('\rho')


% ##        #######   ######  ##     ##  ######  ######## 
% ##       ##     ## ##    ## ##     ## ##    ##    ##    
% ##       ##     ## ##       ##     ## ##          ##    
% ##       ##     ## ##       ##     ##  ######     ##    
% ##       ##     ## ##       ##     ##       ##    ##    
% ##       ##     ## ##    ## ##     ## ##    ##    ##    
% ########  #######   ######   #######   ######     ##    


% load data
load('/local-data/DA-paper/locust/example-data.mat')

% clean up, sub-sample to 1ms
PID = PID1; clear PID1
EAG = EAG1; clear EAG1 

PID = PID(:,1:10:end)';
EAG = EAG(:,1:10:end)';
valve = ODR1(:,1:10:end)';
valve(valve<max(max(valve))/2) = 0;
valve(valve>0) = 1;

% set zero
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:300,i));
	EAG(:,i) = EAG(:,i) - mean(EAG(1:300,i));
end

% filter
PID = bandPass(PID,Inf,30);
EAG = bandPass(EAG,2e3,Inf);
t = 1e-3*(1:length(PID));

[K, EAG_prediction] = extractFilters(PID,EAG,'filter_length',2e3,'filter_offset',500);

% find the gains when the valve turns on
shat = computeSmoothedStimulus(mean(PID,2),500);
[whiff_starts,whiff_ends] = computeOnsOffs(valve(:,1));
rm_this = ((whiff_ends-whiff_starts)<50) | whiff_ends > 29e3 | whiff_starts < 1e3;
whiff_starts(rm_this) = [];
whiff_ends(rm_this) = [];
mean_stim = NaN*whiff_ends;
gain = NaN*whiff_ends;
gain_err =  NaN*whiff_ends;
for i = 1:length(whiff_ends)
	mean_stim(i) = mean(shat(whiff_starts(i):whiff_ends(i)));
	x = mean(EAG_prediction(whiff_starts(i):whiff_ends(i),:),2);
	y = mean(EAG(whiff_starts(i):whiff_ends(i),:),2);
	ff=fit(x(:),y(:),'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_err(i) = diff(temp(:,1))/2;
end
rm_this = gain < 0; 
gain(rm_this) = [];
mean_stim(rm_this) = [];
whiff_starts(rm_this) = [];
whiff_ends(rm_this) = [];


subplot(2,5,5), hold on
plot(mean_stim,gain,'k+');
set(gca,'XScale','log','YScale','log')
xlabel('\mu_{Stimulus} in preceding 500ms (V)')
ylabel('EAG Gain (mV/V)')

% compute rho for various history lengths
tau_h = round(logspace(1,4,50));
rho = computeRhoForGainTimescales(gain,whiff_starts,mean(PID,2),tau_h);
subplot(2,5,10), hold on 
plot(tau_h,rho,'k+')
set(gca,'XScale','log','XTick',[10 1e2 1e3 1e4])
xlabel('History Length (ms)')
ylabel('\rho')

%  ######     ###    ########  ##        #######  ######## ########    ###    
% ##    ##   ## ##   ##     ## ##       ##     ##    ##       ##      ## ##   
% ##        ##   ##  ##     ## ##       ##     ##    ##       ##     ##   ##  
% ##       ##     ## ########  ##       ##     ##    ##       ##    ##     ## 
% ##       ######### ##   ##   ##       ##     ##    ##       ##    ######### 
% ##    ## ##     ## ##    ##  ##       ##     ##    ##       ##    ##     ## 
%  ######  ##     ## ##     ## ########  #######     ##       ##    ##     ## 

% ########     ###    ########    ###    
% ##     ##   ## ##      ##      ## ##   
% ##     ##  ##   ##     ##     ##   ##  
% ##     ## ##     ##    ##    ##     ## 
% ##     ## #########    ##    ######### 
% ##     ## ##     ##    ##    ##     ## 
% ########  ##     ##    ##    ##     ## 




prettyFig('fs',12)

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;



