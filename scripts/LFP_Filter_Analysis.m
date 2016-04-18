

pHeader;


%% LFP_Filter_Analysis
% In this document, we analyse various ways of extracting the filter from the LFP. The problem we are trying to grapple with is that the LFP has some slow drift, and we don't know what the absolute value of the LFP means. There are three ways we could, in principle, extract the filter for the stimulus to LFP transformation:
% 
% # Directly back out filter from stimulus and LFP data
% # Back out filter from LFP and stimulus, but filter the LFP to remove slow fluctuations
% # Back out filter from stimulus and the derivative of the LFP, and integrate that filter
%

%%
% In this document we compare these three methods for all the data we have.

%% Gaussian Stimulus Data
% First, let's see how this comparison plays out for the Gaussian stimulus data. 


[PID, LFP, fA, paradigm,~, ~, AllControlParadigms] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);


% sort the paradigms sensibly
sort_value = [];
for i = length(AllControlParadigms):-1:1
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);

AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)) |  isnan(sum(LFP)));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];

% differentiate the LFP
dLFP = LFP;
for i = 1:length(paradigm)
	% first high pass them to remove spikes
	dLFP(:,i) = bandPass(LFP(:,i),Inf,10);
	dLFP(:,i) = -1e4*filtfilt(ones(10,1),10,[0; diff(dLFP(:,i))]);
end

% filter the LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = filtered_LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,filtered_LFP(:,i));
	filtered_LFP(:,i) = filtered_LFP(:,i)*10; % to get the units right, now in mV
end

% extract filters in all cases
a = 35e3; z = 55e3;
[K1,K1p] = extractFilters(PID,LFP,'use_cache',true,'a',a,'z',z);
[K2,K2p] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);
[K3,K3p] = extractFilters(PID,dLFP,'use_cache',true,'a',a,'z',z);
ft = 1e-3*(1:length(K1)) - .1;

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
subplot(2,2,1), hold on
plot(ft,K1(:,paradigm==1)/norm(K1(:,paradigm==1)),'k')
title('Stimulus \rightarrow LFP')
xlabel('Filter lag (s)')

subplot(2,2,2), hold on
plot(ft,K2(:,paradigm==1)/norm(K2(:,paradigm==1)),'r')
title('Stimulus \rightarrow filtered LFP')
xlabel('Filter lag (s)')

subplot(2,2,3), hold on
temp = -cumsum(K3(:,paradigm==1));
temp = temp/norm(temp);
plot(ft,temp,'b')
title('\int{(Stimulus \rightarrow  dLFP)}')
xlabel('Filter lag (s)')

subplot(2,2,4), hold on
plot(ft,mean(K1(:,paradigm==1)/norm(K1(:,paradigm==1)),2),'k')
plot(ft,mean(K2(:,paradigm==1)/norm(K2(:,paradigm==1)),2),'r')
plot(ft,mean(temp,2),'b')
xlabel('Filter lag (s)')
title('Comparison of the three methods')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% Naturalistic Stimulus Data
% Now we repeat this analysis for the naturalistic stimulus dataset. 

[PID, LFP, fA, paradigm, orn, ~, AllControlParadigms] = consolidateData('/local-data/DA-paper/natural-flickering/with-lfp/ab3/',1);


% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)) |  isnan(sum(LFP)));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

% differentiate the LFP
dLFP = LFP;
for i = 1:length(paradigm)
	% first high pass them to remove spikes
	dLFP(:,i) = bandPass(LFP(:,i),Inf,30);
	dLFP(:,i) = -1e4*filtfilt(ones(10,1),10,[0; diff(dLFP(:,i))]);
end

% filter the LFP
filtered_LFP = LFP;
for i = 1:width(LFP)
	filtered_LFP(:,i) = filtered_LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,filtered_LFP(:,i));
	filtered_LFP(:,i) = filtered_LFP(:,i)*10; % to get the units right, now in mV
end

% extract filters in all cases
a = 15e3; z = 55e3;
[K1,K1p] = extractFilters(PID,LFP,'use_cache',true,'a',a,'z',z);
[K2,K2p] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);
[K3,K3p] = extractFilters(PID,dLFP,'use_cache',true,'a',a,'z',z);
ft = 1e-3*(1:length(K1)) - .1;

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
subplot(2,2,1), hold on
plot(ft,K1(:,orn<4)/norm(K1(:,orn<4)),'k')
title('Stimulus \rightarrow LFP')
xlabel('Filter lag (s)')

subplot(2,2,2), hold on
plot(ft,K2(:,orn<4)/norm(K2(:,orn<4)),'r')
title('Stimulus \rightarrow filtered LFP')
xlabel('Filter lag (s)')

subplot(2,2,3), hold on
temp = -cumsum(K3(:,orn<4));
temp = temp/norm(temp);
plot(ft,temp,'b')
title('\int{(Stimulus \rightarrow  dLFP)}')
xlabel('Filter lag (s)')

subplot(2,2,4), hold on
plot(ft,mean(K1(:,orn<4)/norm(K1(:,orn<4)),2),'k')
plot(ft,mean(K2(:,orn<4)/norm(K2(:,orn<4)),2),'r')
plot(ft,mean(temp,2),'b')
xlabel('Filter lag (s)')
title('Comparison of the three methods')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% 
% Now, we compare the predictions of the LFP using the various filters. In the following figure, we do this for one trial as an illustration. 

% remove baseline from all the LFP
for i = 1:width(LFP)
	LFP(:,i) = LFP(:,i) - mean(LFP(1:5e3,i));
end

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
subplot(2,2,1), hold on
plot(-K1p(a:z,1),-LFP(a:z,1),'k')
title('Stimulus \rightarrow LFP')
xlabel('K_1 \otimes s(t)')
ylabel('LFP (mV)')
legend(['r^2 = ' oval(rsquare(-K1p(a:z,1),-LFP(a:z,1)))],'Location','southeast')

subplot(2,2,2), hold on
plot(-K2p(a:z,1),-LFP(a:z,1),'r')
title('Stimulus \rightarrow filtered LFP')
xlabel('K_2 \otimes s(t)')
ylabel('LFP (mV)')
legend(['r^2 = ' oval(rsquare(-K2p(a:z,1),-LFP(a:z,1)))],'Location','southeast')

subplot(2,2,3), hold on
plot(-K3p(a:z,1),-dLFP(a:z,1),'b')
title('Stimulus \rightarrow dLFP/dt')
xlabel('K_3 \otimes s(t)')
ylabel('dLFP/dt (mV/s)')
legend(['r^2 = ' oval(rsquare(-K3p(a:z,1),-dLFP(a:z,1)))],'Location','southeast')

subplot(2,2,4), hold on
plot(-cumsum(K3p(a:z,1)),-LFP(a:z,1),'b')
title('Stimulus \rightarrow dLFP/dt')
xlabel('\int K_3 \otimes s(t)')
ylabel('LFP (mV)')
legend(['r^2 = ' oval(rsquare(-cumsum(K3p(a:z,1)),-LFP(a:z,1)))],'Location','southeast')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% Integrating the predicted dLFP using the derivative-taking filter is really bad because the derivative-taking filter isn't perfectly derivative taking. So there is a constant trend in the data (since there is a constant offset in the derivative).


%% Version Info
%
pFooter;



