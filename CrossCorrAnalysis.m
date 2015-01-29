% CrossCorrAnalysis.m
% 
% 
% created by Srinivas Gorur-Shandilya at 9:49 , 28 January 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%%
% The purpose of this document is to understand how the cross correlation function relates to the filter. 

calling_func = dbstack;
being_published = 0;
if ~isempty(calling_func)
	if find(strcmp('publish',{calling_func.name}))
		being_published = 1;
	end
end


% load data and prep
load('/local-data/DA-paper/large-variance-flicker/2015_01_22_CS_F1_ab3_3_EtAc.mat')
time = 1e-4*(1:length(data(6).PID));
[fA,tA] = spiketimes2f(spikes(6).A,time);
fA = mean2(fA);
PID = interp1(time,mean2(data(6).PID),tA);


% make up a filter
p.  tau1= 3*5.0122;
p.   K_n= 4.4572;
p.  tau2= 3*20.3750;
p.     A= 41.3682;
p.     n= 4.4297;
p.    Kd= 20.0591;
p.offset= 20.0278;
p.   K_A= 0.3491;

t = 1:1e3;
K = filter_gamma2(p.tau1,p.K_n,p.tau2,p.K_A,t);
t = t*mean(diff(tA));
fp = filter(K,1,PID);

% trivial scaling
fp = fp + 65.5413;
fp = fp*0.1492;

[xc,tc]=xcorr(fp-mean(fp),PID-mean(PID),'unbiased');
xc = xc/max(xc);
tc = tc*mean(diff(tA));

Kx=(xc(tc>-0.1 & tc < 1));
filtertime = (tc(tc>-0.1 & tc < 1));

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(t,K,'k')
plot(filtertime,Kx,'r')





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