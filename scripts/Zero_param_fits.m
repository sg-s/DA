% Zero_param_fits.m
% 
% created by Srinivas Gorur-Shandilya at 1:52 , 10 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.



pHeader;

%% Zero-parameter fits to data
% We attempt to directly measure the parameters of a "zero-parameter" model that includes gain control terms that are sensitive to the mean and the contrast of the stimulus. The model response is given by
% 
% $$ \hat{R}(t)=f\left(g(t)K_{r}\otimes s(t)\right) $$
%
% where 
% 
% $$ g(t)=\frac{1}{1+\beta_{\mu}K_{\mu}\otimes s(t)} $$
%
% is the time-dependent gain of the transduction currents and the nonlinear function is a Hill function. The exponent of the Hill function is also controlled by the stimulus as follows:
%
% $$ n(t)=\frac{n_{0}}{1+\beta_{\sigma}\left|K_{\sigma}\otimes s'(t)\right|} $$ 
%  

%% Fitting a DA Model to the LFP data
% In this section, we fit a modified DA model to the LFP responses to Gaussians in with increasing means. Since we observed that this data shows a Weber-Fechner like gain change, and since we know the DA model can do this, it is reasonable to fit a DA model here. In this section, the DA model we fit is modified so that we remove the mean (since we do the same, essentially for the LFP). 


[PID, LFP, fA, paradigm,~, ~, AllControlParadigms] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);


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

% band pass all the LFP
try 
	load('/local-data/DA-paper/LFP-MSG/september/filtered_LFP.mat','filtered_LFP')
catch
	filtered_LFP = LFP;
	for i = 1:width(LFP)
		filtered_LFP(:,i) = 10*bandPass(LFP(:,i),1e4,Inf);
	end
end

% define limits on data
a = 25e3; z = 45e3;

% fit a DA model here
clear p_DA_LFP d 
for i = 1:max(paradigm)
	d(i).stimulus = nanmean(PID(a:z,paradigm==i),2);
	d(i).response = -nanmean(filtered_LFP(a:z,paradigm==i),2);
	d(i).response(1:1e3) = NaN;
end

p_DA_LFP.   s0 = -0.4541;
p_DA_LFP.  n_z = 2;
p_DA_LFP.tau_z = 46.7500;
p_DA_LFP.  n_y = 1;
p_DA_LFP.tau_y = 76.5000;
p_DA_LFP.    C = 0;
p_DA_LFP.    A = 6.1211;
p_DA_LFP.    B = 1.8301;


%%
% The parameters of the best-fit DA model to this data are:

disp(p_DA_LFP)

%% Fitting a DA Model to the firing data
% For completeness, we also fit a DA model to the firing data in the same experiment (the Weber-Fechner experiment). 

% fit a DA model here
clear p_DA_fA d 
for i = 1:max(paradigm)
	d(i).stimulus = nanmean(PID(a:z,paradigm==i),2);
	d(i).response = nanmean(fA(a:z,paradigm==i),2);
	d(i).response(1:1e3) = NaN;
end

p_DA_fA. s0 = -0.1320;
p_DA_fA.  n_z = 2;
p_DA_fA.tau_z = 135.1249;
p_DA_fA.  n_y = 2;
p_DA_fA.tau_y = 22.2783;
p_DA_fA.    C = 0;
p_DA_fA.    A = 365.4908;
p_DA_fA.    B = 10.4859;

%%
% The parameters of the DA model fit to this data are:
disp(p_DA_fA)

%% Directly estimating extent of gain control from the data
% In this section we attempt to directly measure the gain-control parameters. Since this is an over-determined problem (we have dozens of trials to fit 2 parameters), we find the solution that minimizes the L-2 norm of the error.  

% extract filters and find gain
a = 35e3; z = 55e3;
mean_stim = nanmean(PID(a:z,:));
[K,fA_pred,fA_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);

% compute filter only from the lowest dose
K = K(:,paradigm==1);
K = nanmean(K(:,[2 5]),2);
% use this to make projections everywhere
filtertime = 1e-3*(1:(length(K))) - .1;
fp = NaN*fA;
for i = 1:width(PID)
	fp(:,i) = convolve(1e-3*(1:length(PID)),PID(:,i),K,filtertime);
end

G = 1./fA_gain;
M = ones(length(G),2);
M(:,2) = mean_stim(:);
% remove NaNs
rm_this = isnan(G) | isnan(M(:,2));
G(rm_this) = [];
M(rm_this,:) = [];

X = M\G;
A = 1/X(1); 
B = X(2)*A;

%%
% The A and B parameters estimated through this direct method are:

A,B

%%
% In the following figure, we plot the linear projections corrected by this zero-parameter gain scaling term. Here, we neglect the dynamical nature of the gain control, and assume that the gain simply depends on the mean stimulus over the entire trial. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

subplot(1,2,1), hold on
ss = 100;
c = parula(max(paradigm)+1);

for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(fA_pred(a:z,paradigm == i),2);
	x = x - nanmean(x);
	s = nanmean(PID(a:z,paradigm==i),2);
	x = x + nanmean(s);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('Projected Stimulus (V)')
ylabel('Firing Rate (Hz)')

XG = PID;
for i = 1:width(fA_pred)
	s = nanmean(PID(a:z,i),2);
	x = fA_pred(:,i);
	x = x - nanmean(x);
	x = x + nanmean(s);
	XG(:,i) = (A*x)./(1 + B*mean_stim(i));
end

subplot(1,2,2), hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = nanmean(XG(a:z,paradigm == i),2);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('Gain corrected projection')
suptitle('Gain correction by mean stimulus')
prettyFig('plw',1.3,'lw',1.5,'fs',14,'FixLogX',true,'FixLogY',false)

if being_published
	snapnow
	delete(gcf)
end

%%
% So it seems to be doing pretty well. 

%% Correcting for mean gain control in responses to natural stimuli
% In this section, we attempt to use the parameters we measured in the Weber-Fechner experiment to correct the linear prediction of the response to naturalistic stimuli. 

% load the data
if ~exist('od','var')
	p = '/local-data/DA-paper/natural-flickering/with-lfp/ab3/';
	od = raw2ORNData(p,'filter_LFP',true);

	% specify interesting range
	uts = zeros(length(od(1).stimulus),1);
	uts(10e3:end-5e3) = true;
	for i = 1:length(od)
		od(i).use_this_segment = uts;
	end

	% back out all filters
	od(2) = backOutFilters(od(2));
end

%%
% In the following figure, we plot the response to naturalistic stimuli vs.:
% 
% # the projected stimulus
% # the projected stimulus corrected by dynamic gain measured by fitting a DA model to the LFP in the Weber-Fechner experiment 
% # the projected stimulus corrected by dynamic gain measured by fitting a DA model to the firing rate in the Weber-Fechner experiment
% # the projected stimulus corrected by dynamic gain measured by directly estimating gain parameters in the Weber-Fechner experiment, and using a reasonable timescale
% # the projected stimulus corrected by a dynamic gain term fit to the naturalistic stimulus data. This acts as a positive control. 
% 


example_data = od(2);
orn_data = ORNData;
orn_data.stimulus = nanmean(example_data.stimulus,2);
orn_data.stimulus = orn_data.stimulus - nanmean(orn_data.stimulus(1:5e3));
orn_data.firing_rate = nanmean(example_data.firing_rate,2);
orn_data = backOutFilters(orn_data);
stim_on = false(length(example_data.stimulus),1);
stim_on(10e3:end-5e3) = true; 
orn_data.use_this_segment = stim_on;
orn_data.firing_projected(orn_data.firing_projected<0) = 0;

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,3,1); hold on
x = orn_data.firing_projected;
plot(x(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(x(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')

% DA model fit to LFP Weber data
% correct for some trivial scaling parameters
p_DA_LFP.s0 = -.29;
R = DAModelv2(orn_data.stimulus,p_DA_LFP);
subplot(2,3,2); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Corr. Proj. Stimulus (V)')
title(['Weber-LFP-DA,\tau = ' oval(p_DA_LFP.n_z*p_DA_LFP.tau_z) ,'ms'])

% DA Model fit to Weber firing data
% correct for some trivial scaling parameters
p_DA_fA.s0 = -.06;
R = DAModelv2(orn_data.stimulus,p_DA_fA);
subplot(2,3,3); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Corr. Proj. Stimulus (V)')
title(['Weber-Firing-DA,\tau = ' oval(p_DA_fA.n_z*p_DA_fA.tau_z) ,'ms'])

clear p
p.A = A; p.B = B;
p.n = p_DA_fA.n_z; p.tau = p_DA_fA.tau_z;
S = [orn_data.firing_projected-nanmean(orn_data.firing_projected(1:5e3)), orn_data.stimulus-nanmean(orn_data.stimulus(1:5e3))];
R = adaptiveGainModel(S,p);
subplot(2,3,4); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Corr. Proj. Stimulus (V)')
title(['Weber-directly measured,\tau = ' oval(p_DA_fA.n_z*p_DA_fA.tau_z) ,'ms'])


% best fit DA model
clear p
p.   s0 = 0;
p.  n_z = 1.6015;
p.tau_z = 33.5624;
p.  n_y = 6.9375;
p.tau_y = 6.9189;
p.    C = 0;
p.    A = 450.8111;
p.    B = 5.9000;
[R,~,~,Ky,Kz] = DAModelv2(orn_data.stimulus,p);
subplot(2,3,5); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Corr. Proj. Stimulus (V)')
title(['Best-fit DA,\tau = ' oval(p.n_z*p.tau_z) ,'ms'])

% best fit gain-correction term
clear p d
p.  A = 340.29;
p.  B = 2.24;
p.tau = 28.02;
p.  n = 1;
d.stimulus = [orn_data.firing_projected-min(orn_data.firing_projected) ,orn_data.stimulus];
[R,Kg] = adaptiveGainModel(d.stimulus,p);
subplot(2,3,6); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Corr. Proj. Stimulus (V)')
title(['Best-fit Gain corrected,\tau = ' oval(p_DA_LFP.n_z*p_DA_LFP.tau_z) ,'ms'])

prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end


%% Adding on a contrast-sensitive Hill function
% For each of these models, we fit a contrast-sensitive Hill function as described in other documents to see how they improve the fit. 

clear d p S
S = [0; diff(orn_data.stimulus)];
S_diff = filtfilt(ones(10,1),10,S); clear S
d.stimulus = [orn_data.firing_projected,S_diff];
d.response = orn_data.firing_rate;
d.response(1:1e4) = NaN;

p. n0 = 0.5713;
p.tau = 20.0625;
p.  K = 2.8984;
p.  A = 146.2578;
p.  B = 1.8411e+03;
p.  n = 4.7383;

figure('outerposition',[0 0 1200 700],'PaperUnits','points','PaperSize',[1200 700]); hold on
subplot(2,3,1); hold on
[R,~,n] = contrastLNModel(d.stimulus,p);
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')
title(['Linear Prediction, \tau_{\sigma} = ' oval(p.n*p.tau) 'ms'])

% DA Model LFP Weber
p_DA_LFP.s0 = -.29;
d.stimulus(:,1)  = DAModelv2(orn_data.stimulus,p_DA_LFP);
clear p
p. n0 = 0.5479;
p.tau = 21.0625;
p.  K = 28.8984;
p.  A = 157.5703;
p.  B = 1.9051e+03;
p.  n = 4.9570;
subplot(2,3,2); hold on
R = contrastLNModel(d.stimulus,p);
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')
title(['Weber-LFP-DA, \tau_{\sigma} = ' oval(p.n*p.tau) 'ms'])

% DA Model Weber firing rate
p_DA_fA.s0 = -.06;
d.stimulus(:,1) =  DAModelv2(orn_data.stimulus,p_DA_fA);
clear p
p. n0 = 0.4180;
p.tau = 21.5000;
p.  K = 1.6769e+03;
p.  A = 142.5391;
p.  B = 2.3851e+03;
p.  n = 4.3047;
subplot(2,3,3); hold on
R = contrastLNModel(d.stimulus,p);
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')
title(['Weber-Firing-DA, \tau_{\sigma} = ' oval(p.n*p.tau) 'ms'])


% adaptive gain model
clear p
p.A = A; p.B = B;
p.n = p_DA_fA.n_z; p.tau = p_DA_fA.tau_z;
S = [orn_data.firing_projected-nanmean(orn_data.firing_projected(1:5e3)), orn_data.stimulus-nanmean(orn_data.stimulus(1:5e3))];
d.stimulus(:,1) = adaptiveGainModel(S,p);
clear p
p. n0 = 0.9649;
p.tau = 36.5000;
p.  K = 390.7000;
p.  A = 137.1094;
p.  B = 2.9024e+03;
p.  n = 2.6641;
subplot(2,3,4); hold on
R = contrastLNModel(d.stimulus,p);
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')
title(['Weber-Directly measured, \tau_{\sigma} = ' oval(p.n*p.tau) 'ms'])

% best fit DA model
clear p
p.   s0 = 0;
p.  n_z = 1.6015;
p.tau_z = 33.5624;
p.  n_y = 6.9375;
p.tau_y = 6.9189;
p.    C = 0;
p.    A = 450.8111;
p.    B = 5.9000;
d.stimulus(:,1)  = DAModelv2(orn_data.stimulus,p);
clear p
p.n0 = 1.7656;
p.tau= 18.5000;
p. K = 73.4764;
p. A = 136.0469;
p. B = 1.0855e+03;
p. n = 4.9844;
subplot(2,3,5); hold on
R = contrastLNModel(d.stimulus,p);
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')
title(['Best-fit-DA, \tau_{\sigma} = ' oval(p.n*p.tau) 'ms'])

% best fit gain-correction term
clear p 
p.  A = 340.29;
p.  B = 2.24;
p.tau = 28.02;
p.  n = 1;
d.stimulus(:,1)  = adaptiveGainModel(S,p);
clear p
p. n0 = 0.9219;
p.tau = 17.5000;
p.  K = 201.4764;
p.  A = 134.4531;
p.  B = 1.5975e+03;
p.  n = 4.9844;
subplot(2,3,6); hold on
R = contrastLNModel(d.stimulus,p);
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')
title(['Best-fit-gain-corrected, \tau_{\sigma} = ' oval(p.n*p.tau) 'ms'])

prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


