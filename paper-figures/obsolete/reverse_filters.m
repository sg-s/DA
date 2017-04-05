

pHeader;


clear cdata
cdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
cdata = cleanMSGdata(cdata);

v2struct(cdata)

time = 1e-3*(1:length(PID));

%% Predictive Filters
% In this document, I attempt to back out filters that predict the stimulus given the response. This obviously will not work if the spectral content of the input and the output are very different: for example, white noise inputs can't be predicted from low-passed versions of the white noise, no matter what. To consider a realistic case, I consider the correlated fluctuating stimulus I used to stimualate ab3A ORNs:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
S = PID(:,11);
plot(time,S);

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now, I consider an idealized filter and generate some responses using this. I then flip the time series around and try to back out a filter from the response to the stimulus. 

clear p
p.tau2 = 70;
p.tau1 = 20;
p.   n = 2;
p.   A = 0.3000;
K = filter_gamma2(1:500,p);
R = filter(K,1,S);

a = 15e3;
z = 45e3;
K_rev = fitFilter2Data(flipud(R(a:z)),flipud(S(a:z)),'reg',1,'offset',200,'filter_length',1000);
K_rev(1:100) = [];
filtertime = 1e-3*(1:length(K_rev)) - .1;
sp = convolve(time,flipud(R),K_rev,filtertime);

figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on
plot(K)
xlabel('Lag (ms)')
title('Forward filter')

subplot(1,3,2)
plot(filtertime,K_rev)
xlabel('Anti-lag (s)')
title('Reverse filter')

subplot(1,3,3); hold on
plot(flipud(sp),S,'k.')
title(['r^2 = ' oval(rsquare(flipud(sp),S))])
xlabel('K_{rev} \otimes R')
ylabel('S')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% That seems to work. Now, we apply it to the real data. 

K_rev = NaN(800,length(paradigm));

for i = 1:width(PID)
	if paradigm(i) == 1
		S = PID(:,i);
		R = fA(:,i);
		K_rev(:,i) = fitFilter2Data(flipud(R(a:z)),flipud(S(a:z)),'reg',1,'offset',200,'filter_length',800);
	end
end

K_rev(1:100,:) = [];
filtertime = 1e-3*(1:length(K_rev)) - .1;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(filtertime,nanmean(K_rev,2),'r')
plot(filtertime,nanmean(K2(:,paradigm==1),2),'b')
legend({'K_{rev}','K_{forward}'})
xlabel('Lag (or reverse lag) (s)')
ylabel('Filter amplitude')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% As predicted, the reverse filters are simply the forward filters flipped around. 

%% Filters to higher dimensions of the stimulus. 
% In this section, I try to back out a filter that operates on the square of the derivative of the stimulus, and use it to predict the response. First, let me check if I can reliably estimate the derivative of the stimulus, with some minimal smoothing:

S = PID(:,11);


figure('outerposition',[0 0 1000 699],'PaperUnits','points','PaperSize',[1000 699]); hold on
subplot(2,1,1); hold on
plot(time,S,'k')
plot(time,filtfilt(ones(20,1),20,S),'r')
legend({'Raw Stimulus','Smoothed over 20ms'})
set(gca,'XLim',[20 24])
xlabel('Time (s)')
ylabel('Stimulus (V)')

subplot(2,1,2); hold on
plot(time,[0; diff(S)],'k')
plot(time,[0; diff(filtfilt(ones(20,1),20,S))],'r')
set(gca,'XLim',[20 24])
legend({'Derivative of Raw Stimulus','Derivative of Smoothed stimulus'})
xlabel('Time (s)')
ylabel('dS/dt')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% OK, that works fine. Now, we square this:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ds2 = ([0; diff(filtfilt(ones(20,1),20,S))]).^2;
plot(time,ds2,'r')
xlabel('Time (s)')
ylabel('(ds/dt)^2')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% now we back out filters from this to the response and compare how well it does against the normal filter. 

a = 25e3;
z = 55e3;
K = fitFilter2Data(ds2(a:z),fA(a:z,11),'reg',1,'offset',200,'filter_length',800);
K(1:100) = [];

figure('outerposition',[0 0 601 602],'PaperUnits','points','PaperSize',[601 602]); hold on
subplot(2,2,1); hold on
plot(filtertime,K2(:,11),'k')
title('Stimulus \rightarrow Response')
xlabel('Lag (s)')

subplot(2,2,3); hold on
plot(fA_pred(a:z,11),fA(a:z,11),'k.')
title(['r^2 = ' oval(rsquare(fA_pred(a:z,11),fA(a:z,11)))])
xlabel('K \otimes S')
ylabel('R')

subplot(2,2,2); hold on
plot(filtertime,K,'r')
title('(ds/dt)^2 \rightarrow Response')
xlabel('Lag (s)')

subplot(2,2,4); hold on
fp = convolve(time,ds2,K,filtertime);
plot(fp(a:z),fA(a:z,11),'r.')
title(['r^2 = ' oval(rsquare(fp(a:z),fA(a:z,11)))])
xlabel('K_{diff} \otimes (ds/dt)^2')
ylabel('R')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% 
% Now, we back out filters both to the filter and the squared derivative, at the same time. In the following figure, $K_{1}$ and $K_{2}$ refers to the filters backed out from the stimulus and the squared derivatives simultaneously, and $K_{linear}$ is the normal linear filter we back out from just the stimulus. We also compare the quality of the prediction from these two ways. 

if ~exist('K_dual','var')
	K_dual = NaN(700,length(paradigm),2);

	for i = 1:width(PID)

		if paradigm(i) == 1
			S = PID(:,i);
			ds2 = ([0; diff(filtfilt(ones(20,1),20,S))]).^2;
			R = fA(:,i);
			temp = data2Filter([S ds2],R,'reg',1,'offset',200,'filter_length',900,'left_trim',100,'right_trim',100);
			filtertime = temp.time;
			K_dual(:,i,1) = temp.K(:,1);
			K_dual(:,i,2) = temp.K(:,2);
		end
	end
end

% make predictions from linearly combining these two 
dual_pred = NaN*fA_pred;

for i = 1:length(paradigm)
	if paradigm(i) == 1
		S = PID(:,i); S = S - nanmean(S); S = S/nanstd(S);
		ds2 = ([0; diff(filtfilt(ones(20,1),20,S))]).^2;
		ds2 = ds2 - nanmean(ds2); ds2 = ds2/nanstd(ds2);
		pred1 = convolve(time,S,K_dual(:,i,1),filtertime);
		pred2 = convolve(time,ds2,K_dual(:,i,2),filtertime);
		dual_pred(:,i) = pred1 + pred2;
	end
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(filtertime,nanmean(K_dual(:,:,1),2),'b')
plot(filtertime,nanmean(K_dual(:,:,2),2),'r')
plot(filtertime,nanmean(K2(:,paradigm==1),2),'k')
legend({'K_{1}','K_{2}','K_{linear}'})
xlabel('Lag (s)')
ylabel('Filter amplitude')

subplot(1,2,2); hold on
for i = 1:length(paradigm)
	if paradigm(i) == 1
		x = rsquare(fA_pred(a:z,i),fA(a:z,i));
		y = rsquare(dual_pred(a:z,i),fA(a:z,i));
		plot(x,y,'r+');
	end
end
plot([0.8 1],[0.8 1],'k--')
xlabel('r^2 Linear Model')
ylabel('r^2 Dual filters')


prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Natural Stimuli
% Now I try to do the same for the naturalistic stimulus data
 

load(getPath(dataManager,'5c7dacc5b42ff0eebb980d80fec120c3'),'data','spikes')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;

% A spikes --> firing rate
fA = spiketimes2f(all_spikes,time);

tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 

S = PID - mean(mean(PID(1:5e3,:))); S = mean(S,2);
R = mean(fA,2);
time = 1e-3*(1:length(PID));

% back out linear filter
[K, filtertime_full] = fitFilter2Data(S,R,'reg',1,'filter_length',1999,'offset',500);
filtertime_full = filtertime_full*mean(diff(tA));
filtertime = 1e-3*(-200:900);
K = interp1(filtertime_full,K,filtertime);

% now back out dual filter
ds2 = ([0; diff(filtfilt(ones(20,1),20,S))]).^2;
temp = data2Filter([S ds2],R,'reg',1,'offset',300,'filter_length',1301,'left_trim',100,'right_trim',100);
clear K_dual
K_dual(:,1) = temp.K(:,1);
K_dual(:,2) = temp.K(:,2);


% make linear and dual predictions
fp = convolve(time,S,K,filtertime);

stim = S - nanmean(S); stim = stim/nanstd(stim);
ds2 = ([0; diff(filtfilt(ones(20,1),20,stim))]).^2;
ds2 = ds2 - nanmean(ds2); ds2 = ds2/nanstd(ds2);
pred1 = convolve(time,stim,K_dual(:,1),filtertime);
pred2 = convolve(time,ds2,K_dual(:,2),filtertime);
dual_pred = pred1 + pred2;

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
plot(filtertime,K,'k')
plot(filtertime,K_dual(:,1),'b')
plot(filtertime,K_dual(:,2),'r')
legend({'K_{linear}','K_1','K_2'})

subplot(1,3,2); hold on
plot(fp,R,'k.')
xlabel('Linear Prediction')
ylabel('Response (Hz)')

subplot(1,3,3); hold on
plot(dual_pred,R,'k.')
xlabel('Dual Prediction')
ylabel('Response (Hz)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;

