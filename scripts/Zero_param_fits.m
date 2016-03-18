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

%% 
% The ultimate goal here is to use this model to account for ORN responses to natural stimuli. We will do so in the following stages:
% 
% # Fitting a DA Model to the data. The term inside the Hill function in the "zero-parameter" model is essentially a DA model. 
% # Correcting the linear prediction of the naturalistic stimulus by a divisive gain term. This is identical to the previous case, but we directly measure the response filter from the data, instead of fitting it. 
% 

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

%% 1. Fitting a DA Model to the data
% First, we fit a DA Model to the data. The following figure shows the ORN response plotted vs. the best-fit DA model, and the resultant goodness of fit. The point of this section is to make sure that a DA model fits the data well, and outperforms LN models, etc. 

% do the analysis of fast gain control
example_data = od(2);
orn_data = ORNData;
orn_data.stimulus = nanmean(example_data.stimulus,2);
orn_data.stimulus = orn_data.stimulus - nanmean(orn_data.stimulus(1:5e3));
orn_data.firing_rate = nanmean(example_data.firing_rate,2);
orn_data = backOutFilters(orn_data);
stim_on = false(length(example_data.stimulus),1);
stim_on(10e3:end-5e3) = true; 
orn_data.use_this_segment = stim_on;

% also fit a NL
orn_data_LN = fitNL(orn_data);

clear p
p.   s0 = 0;
p.  n_z = 2.6562;
p.tau_z = 94.2499;
p.  n_y = 10.5938;
p.tau_y = 4.0127;
p.    C = 0.6004;
p.    A = 712.5611;
p.    B = 11.8219;
[R,~,~,K_DA,K_DA2] = DAModelv2(orn_data.stimulus,p);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
x = orn_data.firing_projected;
[~,excursions] = plotExcursions(orn_data,gca);
cla
plot(x(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')

% make two time vectors, one defining when the stimulus is on, and one just for the whiffs
whiff_times = false(length(orn_data.stimulus),1);
for i = 1:length(excursions.ons)
	whiff_times(excursions.ons(i):excursions.offs(i)) = true;
end
legend(['r^2 = ' oval(rsquare(x(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')

% fit a NL just to the excursions
orn_data.use_this_segment = whiff_times;
orn_data_LN = fitNL(orn_data,'firing_rate');
orn_data_LN.use_this_segment = stim_on;

subplot(1,3,2); hold on
x = orn_data_LN.firing_projected;
plot(x(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(x(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('LN Model Prediction (Hz)')

subplot(1,3,3); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('DA Model Prediction (Hz)')

prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% It looks like it does pretty OK. What do the DA model filters look like?

figure('outerposition',[0 0 600 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(K_DA)
plot(K_DA2)
xlabel('Filter Lag (ms)')
legend('K_{y}','K_{z}')
set(gca,'XLim',[0 600],'YLim',[0 0.03])
prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% The parameters of this DA model are:

disp(p)


%%
% Can we fit a DA model keeping C = 0? This would be equivalent to the inner term of the "zero parameter" model. 

clear p
p.   s0 = 0;
p.  n_z = 1.6015;
p.tau_z = 33.5624;
p.  n_y = 6.9375;
p.tau_y = 6.9189;
p.    C = 0;
p.    A = 450.8111;
p.    B = 5.9000;
[R,~,~,K_DA,K_DA2] = DAModelv2(orn_data.stimulus,p);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
plot(K_DA)
plot(K_DA2)
xlabel('Filter Lag (ms)')
legend('K_{y}','K_{z}')
set(gca,'XLim',[0 400],'YLim',[0 0.03])


subplot(1,3,2); hold on
x = orn_data_LN.firing_projected;
plot(x(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(x(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('LN Model Prediction (Hz)')

subplot(1,3,3); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('DA Model Prediction (Hz)')

prettyFig('fs',18)


if being_published	
	snapnow	
	delete(gcf)
end

%% 2. Adding a gain-correction term to the linear projection
% Now we attempt to add a gain-correcting term to the linear projection. The motivation here is that we can attempt to directly measure the response filter in the DA model (in the numerator) from the data. We can then add on a mean-sensitive gain term (like the denominator of the DA model) to the linear projection and this should allow us to rebuild the DA model piece by piece. 

clear p d
p.  A = 340.29;
p.  B = 2.24;
p.tau = 28.02;
p.  n = 1;
d.stimulus = [orn_data.firing_projected-min(orn_data.firing_projected) ,orn_data.stimulus];
R = adaptiveGainModel(d.stimulus,p);

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
x = orn_data.firing_projected;
plot(x(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(x(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')

subplot(1,3,2); hold on
x = orn_data_LN.firing_projected;
plot(x(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(x(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('LN Model Prediction (Hz)')

subplot(1,3,3); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Gain-corrected Prediction (Hz)')

prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end


%% 
% The gain-corrected model is barely at the level of the LN model. Clearly, it is not working well. This is clearer when we realize that the best fit gain filter has a timescale of 30ms. Why is this so bad? In theory, this could have been as good as the DA model. In the following figure, we compare the DA response filter to the computed response filter:

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
t = 1e-3*(1:length(K_DA));
plot(t,K_DA/norm(K_DA),'r')
K = orn_data.K_firing(100:end);
t = 1e-3*(1:length(K));
plot(t,K/norm(K),'k')
set(gca,'XLim',[0 .5])
xlabel('Filter Lag (s)')
legend({'DA response filter','Linear filter'})
prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end


%%
% It is quite possible that our filter extracted from this data is pretty bad, since the stimulus is so non-Gaussian. Instead, we attempt to parametrically fit a filter to the data.

clear p
p.tau1 = 25.5000;
p.tau2 = 1;
p.   A = 0;
p.   n = 38.6289;
Kp = filter_gamma2(1:500,p);

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
t = 1e-3*(1:length(K_DA));
plot(t,K_DA/norm(K_DA),'r')
K = orn_data.K_firing(100:end);
t = 1e-3*(1:length(K));
plot(t,K/norm(K),'k')
t = 1e-3*(1:length(Kp));
Kp = filter_gamma2(1:500,p);
plot(t,Kp/norm(Kp),'b')
set(gca,'XLim',[0 .5])
xlabel('Filter Lag (s)')
legend({'DA response filter','Linear filter','Parametric Filter'})
prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% How well does this parametric filter work? In the following figure, we compare the data to parametric filter, the LN model (using a parametric filter) and a gain-corrected linear projection (using the parametric filter). 

p_data = ORNData;
p_data.stimulus = nanmean(example_data.stimulus,2);
p_data.stimulus = p_data.stimulus - nanmean(p_data.stimulus(1:5e3));
p_data.firing_rate = nanmean(example_data.firing_rate,2);
t = 1e-3*(1:length(Kp));
p_data.filtertime_firing = t;
p_data.K_firing = Kp;

% fit a NL just to the excursions
p_data.use_this_segment = whiff_times;
p_data_LN = fitNL(p_data,'firing_rate');
p_data_LN.use_this_segment = stim_on;

figure('outerposition',[0 0 1500 500],'PaperUnits','points','PaperSize',[1500 500]); hold on
subplot(1,3,1); hold on
x = p_data.firing_projected;
plot(x(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(x(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Proj. Stimulus (V)')


subplot(1,3,2); hold on
x = p_data_LN.firing_projected;
plot(x(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(x(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('LN Model Prediction (Hz)')

clear d
d.stimulus = [p_data.firing_projected-min(p_data.firing_projected) ,p_data.stimulus];
d.response = orn_data.firing_rate;
d.response(~whiff_times) = NaN;

clear p
p.  A = 453.3535;
p.  B = 4.8174;
p.tau = 12.7794;
p.  n = 1.1839;
R = adaptiveGainModel(d.stimulus,p);

subplot(1,3,3); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('Gain-corrected Prediction (Hz)')

prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% There is clearly a problem. When we fit both filters together (as we do when we fit a DA model), we get a prediction that is really good. When we try to fit one filter after the other, we get a combination that is much worse than if we fit both filters. To verify that this is true, we set the response filter to be the response filter in the DA model and verify that we can recover the "gain filter" of the DA model: 

clear p_data
p_data = ORNData;
p_data.stimulus = nanmean(example_data.stimulus,2);
p_data.stimulus = p_data.stimulus - nanmean(p_data.stimulus(1:5e3));
p_data.firing_rate = nanmean(example_data.firing_rate,2);
t = 1e-3*(1:length(K_DA));
p_data.filtertime_firing = t;
p_data.K_firing = K_DA;



%% Wrapping a DA Model with a contrast-sensitive term
% In light of these complications, we try a different tack. The DA model is still the model that explains most of the data in this case. In our "zero-parameter" model, we essentially have a DA model wrapped by a Hill function whose steepness changes with the contrast. In this section, we ask if we can directly wrap the DA model prediction with the Hill function, with parameters measured from the contrast switching experiment, to see if we can improve on it. 

%%
% First, we show a proof of concept, fitting the contrast-sensitive part of the model alone around a DA model prediction. To be clear, all the parameters of the contrast-sensitive Hill function were allowed to vary. The best fit parameters were:
clear p
p.   s0 = 0;
p.  n_z = 2.6562;
p.tau_z = 94.2499;
p.  n_y = 10.5938;
p.tau_y = 4.0127;
p.    C = 0.6004;
p.    A = 712.5611;
p.    B = 11.8219;
R = DAModelv2(orn_data.stimulus,p);
R = R - min(R);
R = R/max(R);

clear p
p. n0 = 2.4531;
p.tau = 309;
p.  K = 0.4764;
p.  A = 89.9531;
p.  B = 2.4935e+03;
p.  n = 2.9844;
S = filtfilt(ones(10,1),10,[0; diff(p_data.stimulus)]);
d.stimulus = [R, S];
Rc = contrastLNModel(d.stimulus,p);

disp(p)


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('DA Model (norm)')
ylabel('ORN Response (Hz)')

subplot(1,2,2); hold on
plot(Rc(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(Rc(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel(['DA Model with contrast-' char(10) 'sensitive gain  (Hz)'])
ylabel('ORN Response (Hz)')

prettyFig('fs',18)

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Can we improve on the DA model by plugging in the parameters we measured from the contrast-switch experiment? 

clear p
p. n0 = 7.62;
p.tau = 121/2;
p.  K = 0.355;
p.  A = 66.9531;
p.  B = 1.76e+03;
p.  n = 2;
S = filtfilt(ones(10,1),10,[0; diff(p_data.stimulus)]);
d.stimulus = [R, S];
Rc = contrastLNModel(d.stimulus,p);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(R(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(R(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel('DA Model (norm)')
ylabel('ORN Response (Hz)')

subplot(1,2,2); hold on
plot(Rc(1e4:10:end),orn_data.firing_rate(1e4:10:end),'k.')
legend(['r^2 = ' oval(rsquare(Rc(1e4:10:end),orn_data.firing_rate(1e4:10:end)))],'Location','southeast')
xlabel(['DA Model with contrast-' char(10) 'sensitive gain  (Hz)'])
ylabel('ORN Response (Hz)')

prettyFig('fs',18)
if being_published	
	snapnow	
	delete(gcf)
end


%%
% Clearly, something is wrong, perhaps the timescale. 

%% Version Info
%
pFooter;


