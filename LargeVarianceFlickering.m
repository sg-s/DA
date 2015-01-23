% LargeVarianceFlickering.m
% 
% created by Srinivas Gorur-Shandilya at 3:30 , 19 January 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end

% redo = 1; % deliberately unset

%% Response of ORNs to flickering odor stimuli with large variances
% From previous experiments, we see that the response of ORNs to largely varying stimuli that mimics the natural odor plumes is particularly interesting: in that no model we have can precisely account for the data, and in that, from other results, we expect to see a large variation of gain of the ORNs to these stimuli. 

%%
% In this document, we generate odor stimuli flickers over a large range, like the "natural" stimuli, but never goes to zero, so that the neuron should never silence (allowing us to accurately follow its response). 

%%
% The odor stimulus looks like this (different colours are different standard deviations of the exponentiated Gaussians):

load('/local-data/DA-paper/large-variance-flicker/2015_01_22_CS_F1_ab3_3_EtAc.mat')

figure('outerposition',[0 0 1400 1000],'PaperUnits','points','PaperSize',[1400 1000]); hold on
subplot(2,8,9:14), hold on
time = 1e-4*(1:length(data(6).PID));
plot(time,mean2(data(6).PID),'k')
set(gca,'XLim',[10 60])
xlabel('Time (s)')
ylabel('PID (V)')

subplot(2,8,15:16), hold on
r = rsquare(data(6).PID);
imagescnan(r)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('min r^2=',oval(min(min(r)),2)))

subplot(2,8,1:6), hold on
time = 1e-4*(1:length(data(6).MFC200));
plot(time,mean2(data(6).MFC200),'k')
set(gca,'XLim',[10 60])
xlabel('Time (s)')
ylabel('MFC Flow Signal (V)')

subplot(2,8,7:8), hold on
r = rsquare(data(6).MFC200);
imagescnan(r)
caxis([0 1])
colorbar
axis image
axis off
title(strcat('min r^2=',oval(min(min(r)),2)))

PrettyFig;

if being_published
	snapnow
	delete(gcf)
end


[fA,tA] = spiketimes2f(spikes(6).A,time);
fA = mean2(fA);
PID = interp1(time,mean2(data(6).PID),tA);

%%
% The following figure shows the response of 1 ab3A neuron to this stimulus, with a best fit LN Model. 

p.  tau1= 0.2300;
p.   K_n= 1.4219;
p.  tau2= 39.2500;
p.   K_A= -1.8672;
p.     A= 170.3750;
p.     n= 2;
p.    Kd= -39.2500;
p.offset= 14.0098;

fp = pLNModel(PID,p);


figure('outerposition',[0 0 1000 700],'PaperUnits','points','PaperSize',[1000 700]); hold on
subplot(2,1,1), hold on
plot(tA,PID,'k')
xlabel('Time (s)')
ylabel('PID (V)')
subplot(2,1,2), hold on
plot(tA,fA,'k')
l=plot(tA,fp,'r');
r2 = oval(rsquare(fA,fp),2);
legend(l,strcat('r^2=',r2))
ylabel('Firing Rate (Hz)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% The LN Model used here has been parametrised and looks like this:

t = 1:300;
K = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);
t = t*mean(diff(tA));

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1), hold on
plot(t,K,'r')
xlabel('Filter Lag (s)')
subplot(1,2,2), hold on
plot(1:100,hill([p.A p.Kd p.n],1:100),'r')
ylabel('Predicted f (Hz)')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% How is gain controlled in response to this flickering stimulus? In the following plot we compute the instantaneous gain by dividing the response of the neuron by this LN model prediction. 

gain = fA(:)./fp(:);

figure('outerposition',[0 0 1300 700],'PaperUnits','points','PaperSize',[1300 700]); hold on
subplot(2,1,1), hold on
plot(tA,PID,'k')
ylabel('Stimulus')
subplot(2,1,2), hold on
plot([-1 61],[1 1],'k--')
plot(tA,gain,'r')
ylabel('Gain')

PrettyFig;
if being_published
	snapnow
	delete(gcf)
end

%%
% It looks like the gain is changing ~2 fold in some points. In the following figure, we do a more detailed gain analysis, splitting the data according to when the stimulus is high or low in the past and checking the gain in those points (as before). 

% do gain analysis
clear x
x.response = fA; 
x.prediction = fp;
x.stimulus = PID; 
x.time = tA;
x.filter_length = 299;
ph = [];

history_lengths = (3*floor(1000*logspace(-1.5,1,30)/3))/1e3;
example_history_length = 0.135;

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

if redo
	[p_LN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	s=abs(l-h);
	s(p_LN(1,:)>0.05)=NaN;
	[~,loc]=max(s);

	% save it for later
	ehl = history_lengths(loc);
	pb = p_LN;

else
	GainAnalysis4(x,history_lengths,ehl,ph,pb);
end

xlabel(ph(3),'LN Prediction (Hz)')
set(ph(4),'XScale','log')

	PrettyFig;
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
disp(DataHash(strcat(mfilename,'.m'),Opt))

%%
% This file should be in this commit:
[status,m]=unix('git rev-parse HEAD');
if ~status
	disp(m)
end

