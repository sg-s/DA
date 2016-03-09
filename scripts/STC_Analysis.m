% STC_Analysis.m
% 
% created by Srinivas Gorur-Shandilya at 2:24 , 09 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;

%% STC Analysis
% In this document we explore the possibility of using STC-based methods to pull out two filters that explain ORN response data. We use the iSTAC algorithm made available from the Pillow Lab for computing the STC. 

%% 0: Synthetic Data: Simple LN Model
% In this section we generate synthetic data using a simple LN model (with a quadratic nonlinearity) and attempt to pull out two filters from this. Since there is only one filter that generates this data, we don't expect to be able to pull out a 2nd filter. We expect to see that the 1st Eigenvector of STC accounts for most of the variance, with minimal contributions from adding additional terms. 

% use our own filters
nt = 100;
tvec = (-nt+1:0)';

p.tau2 = 20;
p.tau1 = 5;
p.n = 2;
p.A = .3;
filt1 = fliplr(filter_gamma2(1:100,p))';
filt1 = filt1./norm(filt1);  %normalize
filt2 = 0*filt1;

% % Create synthetic data ------------------
slen = 10000;   % Stimulus length (Better convergence w/ longer stimulus)
Stim = randn(slen,1);
Stim = filtfilt(ones(10,1),10,Stim) + .5;
r = sameconv(Stim,filt1);
r = r.^2;

[K1,K2,~,KLcontributed] = backOutSTC(Stim,r);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
plot(tvec,filt1,'k')
title('1st filter')
plot(tvec,K1,'r')
legend({'Actual','Reconstructed'},'Location','northwest')

subplot(1,3,2); hold on
plot(tvec,filt2,'k')
title('2nd filter')
plot(tvec,K2,'r')

subplot(1,3,3); hold on
plot(1:length(KLcontributed),KLcontributed,'ko')
xlabel('# of dimension')
ylabel('KL contributed')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% 1: Synthetic Data: Two orthogonal filters
% In this section we generate synthetic data using two orthogonal filters, and attempt to back them out using STC. Note that since we have a quadratic non-linearity, some the reconstructed filters may be negations of the actual filters. 

clearvars -except being_published tvec nt
p.tau2 = 20;
p.tau1 = 5;
p.n = 2;
p.A = .3;
filt1 = fliplr(filter_gamma2(1:100,p))';
filt1 = filt1./norm(filt1);  %normalize

filt2 = [diff(filt1); 0];  % 2nd filter
filt2 = filt2- filt1*(filt1'*filt2); %orthogonalize to 1st filter
filt2 = filt2./norm(filt2); % normalize

% % Create synthetic data ------------------
slen = 10000;   % Stimulus length (Better convergence w/ longer stimulus)
Stim = randn(slen,1);

DC = [0.75 .5];  % Setting to zero means expected STA is zero
linresp = [sameconv(Stim, filt1)+DC(1) sameconv(Stim,filt2)+DC(2)];  % filter output
r = 10*linresp(:,1).^2 + 8*linresp(:,2).^2; % instantaneous spike rate

[K1,K2,~,KLcontributed] = backOutSTC(Stim,r);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
plot(tvec,filt1,'k')
title('1st filter')
plot(tvec,K1,'r')
legend({'Actual','Reconstructed'},'Location','northwest')

subplot(1,3,2); hold on
plot(tvec,filt2,'k')
title('2nd filter')
plot(tvec,K2,'r')

subplot(1,3,3); hold on
plot(1:length(KLcontributed),KLcontributed,'ko')
xlabel('# of dimension')
ylabel('KL contributed')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% 2: Synthetic Data: Two non-orthogonal filters
% In this section we generate synthetic data using two non-orthogonal filters, and attempt to back them out using STC. 

clearvars -except being_published tvec nt
p.tau2 = 20;
p.tau1 = 5;
p.n = 2;
p.A = .3;
filt1 = fliplr(filter_gamma2(1:100,p))';
filt1 = filt1./norm(filt1);  %normalize

p.tau2 = 20;
p.tau1 = 7;
p.n = 6;
filt2 = fliplr(filter_gamma2(1:100,p))';
filt2 = filt2./norm(filt2);  %normalize

% Create synthetic data ------------------
slen = 10000;   % Stimulus length (Better convergence w/ longer stimulus)
Stim = randn(slen,1);

DC = [0.75 .5];  % Setting to zero means expected STA is zero
linresp = [sameconv(Stim, filt1)+DC(1) sameconv(Stim,filt2)+DC(2)];  % filter output
r = 10*linresp(:,1).^2 + 8*linresp(:,2).^2; % instantaneous spike rate

[K1,K2,~,KLcontributed] = backOutSTC(Stim,r);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
plot(tvec,filt1,'k')
title('1st filter')
plot(tvec,K1,'r')
legend({'Actual','Reconstructed'},'Location','northwest')

subplot(1,3,2); hold on
plot(tvec,filt2,'k')
title('2nd filter')
plot(tvec,K2,'r')

subplot(1,3,3); hold on
plot(1:length(KLcontributed),KLcontributed,'ko')
xlabel('# of dimension')
ylabel('KL contributed')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% 3: Synthetic Data: DA Model
% In this section we generate synthetic data using the DA Model, and attempt to back them out using STC. We expect to see contributions from more than just the 1st Eigenvector as there truly are more than one filter in the DA model. 

clearvars -except being_published tvec nt
p = '/local-data/DA-paper/natural-flickering/with-lfp/ab3/';
PID = consolidateData(p,true);
S = nanmean(PID,2);

clear p
p.   s0 = 7.8242e-04;
p.  n_z = 2;
p.tau_z = 151.1249;
p.  n_y = 2;
p.tau_y = 26.7002;
p.    C = 0.5457;
p.    A = 163.2252;
p.    B = 2.4703;

% Create synthetic data ------------------
slen = 10000;   % Stimulus length (Better convergence w/ longer stimulus)
Stim = randn(slen,1);

[R,~,~,Ky,Kz] = DAModelv2(S,p);

[K1,K2,~,KLcontributed] = backOutSTC(S(1:10:end),R(1:10:end));

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
plot(Ky)
plot(Kz)
title('DA model filters')
xlabel('Filter Lag (ms)')
ylabel('Filter')

subplot(1,3,2); hold on
plot(tvec*10,K1)
plot(tvec*10,K2)
title('STC filter')

subplot(1,3,3); hold on
plot(1:length(KLcontributed),KLcontributed,'ko')
xlabel('# of dimension')
ylabel('KL contributed')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% 4: Real Data: Gaussian Inputs
% In this section we use real data with Gaussian inputs. Note that we only back out one filter, which makes sense because a simple LN model fits this data well. 

clearvars -except being_published tvec nt
[PID, LFP, fA, paradigm, orn, fly, AllControlParadigms, paradigm_hashes] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

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

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)));
PID(:,bad_trials) = [];
LFP(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];
orn(bad_trials) = [];

S = (PID(25e3:10:55e3,paradigm==2)); S = S(:);
R = (fA(25e3:10:55e3,paradigm==2)); R = R(:);

[K1,K2,~,KLcontributed] = backOutSTC(S,R);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
title('1st filter')
plot(tvec*10,K1,'r')
xlabel('Filter Lag (s)')

subplot(1,3,2); hold on
title('2nd filter')
plot(tvec*10,K2,'r')
xlabel('Filter Lag (s)')

subplot(1,3,3); hold on
plot(1:length(KLcontributed),KLcontributed,'ko')
xlabel('# of dimension')
ylabel('KL contributed')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% 5: Real Data: Naturalistic Stimulus
% In this section we use our naturalistic stimulus data, where stimulus is far from Gaussian. Note that here, in contrast to the real data driven by Gaussian noise, we have contributions from more than one mode, suggesting that there is something beyond a simple LN model at play here. This is consistent with our other results where we show that the LN model leaves some variance unexplained.  

clearvars -except being_published tvec nt
p = '/local-data/DA-paper/natural-flickering/with-lfp/ab3/';
[PID, ~, fA,~,orn] = consolidateData(p,true);
R = fA(1:10:end,orn==4); R = R(:);
S = PID(1:10:end,orn==4); S = S(:);

[K1,K2,STA,KLcontributed] = backOutSTC(S,R);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
title('1st filter')
plot(tvec*10,K1,'r')
xlabel('Filter Lag (ms)')

subplot(1,3,2); hold on
title('2nd filter')
plot(tvec*10,K2,'r')
xlabel('Filter Lag (ms)')

subplot(1,3,3); hold on
plot(1:length(KLcontributed),KLcontributed,'ko')
xlabel('# of dimension')
ylabel('KL contributed')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%% Version Info
%
pFooter;


