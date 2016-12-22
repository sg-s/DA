pHeader;


%% One model to fit all the data 
% In this document I attempt to find one model that fits all the data we have so far. First, I'll fit some simplified version of the model to each data set, and study how well it constrains the parameters of this model. Then, finally, I'll attempt to fit one model to all the data, and compare that model's prediction to the prediction from models that are fit to one dataset alone. 

% first get all the data

% get the filter from the Gaussian stimuli 
clear MSGdata
MSGdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
MSGdata = cleanMSGdata(MSGdata);

% get the sparse naturalistic stimuli with LFP
% get this data
clearvars -except being_published MSGdata
load('/local-data/DA-paper/data-for-paper/nat-stim/ab3A_nat_stim.ORNData','-mat')
od(1) = [];
SNSdata.LFP = [];
SNSdata.fA = [];
SNSdata.PID = [];
SNSdata.orn = [];
for i = length(od):-1:1
	SNSdata.LFP = [SNSdata.LFP  od(i).LFP];
	SNSdata.orn = [SNSdata.orn; vectorise(i*ones(od(i).n_trials,1))];
	SNSdata.fA = [SNSdata.fA  od(i).firing_rate];
	SNSdata.PID = [SNSdata.PID  od(i).stimulus];
end


%% Firing rate: Sparse naturalistic stimulus
% First, I fit a NLN model to the firing rate responses to sparse naturalistic stimuli. The following figure shows the best fit model responses, and I also plot the 

clear p
p.   k0 = 0.1386;
p.tau_z = 10;
p.    B = 0;
p.  n_z = 1;
p.    n = 1;
p. tau1 = 24.2813;
p. tau2 = 78.3750;
p.  n_y = 2;
p.    A = 0.5356;
p.    C = 176.2152;

S = mean(SNSdata.PID(:,SNSdata.orn == 1),2);
R = mean(SNSdata.fA(:,SNSdata.orn == 1),2);
Rhat = aNLN2(S,p);
t = 1e-3*(1:length(R));

figure('outerposition',[0 0 1409 801],'PaperUnits','points','PaperSize',[1409 801]); hold on
subplot(2,4,1:4); hold on
plot(t,R,'k')
plot(t,Rhat,'r')
xlabel('Time (s)')
ylabel('ab3A Firing rate (Hz)')
legend({'ab3A','NLN model'})

subplot(2,4,5); hold on
x = logspace(-4,log10(max(S)*10),100);
a = hill([1 p.k0 p.n],x);
plot(x/mean(S),a,'r')
xlabel('S/<S>')
ylabel('a')
set(gca,'XScale','log','XLim',[1e-3 1e2],'XTick',[1e-3  1e-2 1e-1 1 10])

subplot(2,4,6); hold on
t = 1:600;
clear q
q.A = p.A; q.tau1 = p.tau1; q.tau2 = p.tau2; q.n = p.n_y;
K = filter_gamma2(t,q);
plot(t,K,'r')
xlabel('Filter lag (ms)')
ylabel('Filter')

subplot(2,4,7); hold on
plot(Rhat,R,'k.')
xlabel('NLN prediction (Hz)')
ylabel('ab3A firing rate (Hz)')
legend(['r^2 = ' oval(rsquare(R,Rhat))],'Location','northwest')

% vary k_D and n
all_k_D = logspace(log10(mean(S)/10),log10(mean(S)*10),31);
all_n = 1:10;
r2 = NaN(length(all_k_D),length(all_n));
if exist('.cache/sparse_nat_stim_k_D_n.mat','file')
	load('.cache/sparse_nat_stim_k_D_n.mat')
else
	for i = 1:length(all_k_D)
		textbar(i,length(all_k_D))
		for j = 1:length(all_n)
			q = p;
			q.k0 = all_k_D(i);
			q.n = all_n(j);
			Rhat = aNLN2(S,q);
			r2(i,j) = rsquare(Rhat,R);
		end
	end
	save('.cache/sparse_nat_stim_k_D_n.mat','r2')
end


n_labels = {};
for i = 1:2:length(all_n)
	n_labels{i} = oval(all_n(i));
end
k_D_labels = {'1/10','1','10'};

subplot(2,4,8); hold on
imagesc(r2)
ylabel('K_D/<S>')
xlabel('n')
colorbar
caxis([0.5 1])
set(gca,'XTick',[1:2:10],'XTickLabel',n_labels(cellfun(@(x) ~isempty(x),(n_labels))),'XTickLabelRotation',45)
set(gca,'YTick',[1 16 31],'YTickLabel',k_D_labels)
set(gca,'XLim',[.6 10.5],'YLim',[.5 31.5])
title('r^2(data, NLN model)')


prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Firing rate: dense naturalistic stimulus
%

%% Firing rate: mean shifted Gaussians
% 

%% Firing rate: variance gain control
% 

%% Firing rate: All data
% 


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


