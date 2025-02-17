
pHeader;

%%
% OK, so it looks like both the DA model and the NLN model fit the data well. In the following section, I attempt to fit a NLN model to synthetic data generated by the DA model. The idea is that if we consider "reality" to actually have a fast gain control mechanism in it, can a NLN model reproduce it? In other words, can the stimulus we use distinguish between a DA and a NLN model?

% get the first nat. stim
clear ab3 ab2
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

% remove the baseline from the PID, and remember the error
PID_baseline = mean(mean(PID(1:5e3,:)));
PID = PID - PID_baseline;


% fit a DA model to this
clear p
p.   s0 = 0;
p.  n_z = 2;
p.tau_z = 147.3750;
p.  n_y = 2;
p.tau_y = 27.2500;
p.    C = 0.5000;
p.    A = 170.4375;
p.    B = 2.7656;

S = mean(PID,2);
S = S - min(S);

[R,~,~,Ky,Kz] = DAModelv2(S,p);

clear data
data.response = R;
data.stimulus = [S R];

% fit a NLN model to this
clear p
p.n = 1;
p.k_D = .7769;

R_NLN = NLNmodel(data.stimulus,p);
time = 1e-3*(1:length(R));

%%
% In the following figure, I show the filters of the DA model that generate this data, and use this DA model to generate responses using the sparse naturalistic stimulus. I then fit a NL model to this synthetic data, and show how well the model performs. 

figure('outerposition',[0 0 1300 801],'PaperUnits','points','PaperSize',[1300 801]); hold on
subplot(2,2,1:2); hold on
plot(time,R,'k')
plot(time,R_NLN,'r')
legend({'DA model Response',['NLN model fit, r^2 = ' oval(rsquare(R,R_NLN))]},'Location','northwest')
subplot(2,2,3); hold on
plot(Ky)
plot(Kz)
title('DA model filters')
xlabel('Filter lag (ms)')
ylabel('Filter amplitude')

subplot(2,2,4); hold on
plot(R_NLN,R)
xlabel('NLN model prediction')
ylabel('DA model response')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Maybe this is simply because the gain filter is so close to the "response" filter. What if I modify the gain filter to be different? 


% fit a DA model to this
clear p
p.   s0 = 0;
p.  n_z = 2;
p.tau_z = 70.3750;
p.  n_y = 2;
p.tau_y = 20.2500;
p.    C = 0;
p.    A = 170.4375;
p.    B = 2.7656;

S = mean(PID,2);
S = S - min(S);

[R,~,~,Ky,Kz] = DAModelv2(S,p);

clear data
data.response = R;
data.stimulus = [S R];

% fit a NLN model to this
clear p
p.n = 1;
p.k_D = 2.7574;

R_NLN = NLNmodel(data.stimulus,p);
time = 1e-3*(1:length(R));


figure('outerposition',[0 0 1300 801],'PaperUnits','points','PaperSize',[1300 801]); hold on
subplot(2,2,1:2); hold on
plot(time,R,'k')
plot(time,R_NLN,'r')
legend({'DA model Response',['NLN model fit, r^2 = ' oval(rsquare(R,R_NLN))]},'Location','northwest')
subplot(2,2,3); hold on
plot(Ky)
plot(Kz)
title('DA model filters')
xlabel('Filter lag (ms)')
ylabel('Filter amplitude')

subplot(2,2,4); hold on
plot(R_NLN,R)
xlabel('NLN model prediction')
ylabel('DA model response')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% Version Info
%
pFooter;


