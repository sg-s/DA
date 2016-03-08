% 
% 
% created by Srinivas Gorur-Shandilya at 1:59 , 07 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;

% load the data
if ~exist('od','var')
	p = '/local-data/DA-paper/natural-flickering/with-lfp/ab3/';
	od = raw2ORNData(p);
	% ignore all the LFP for now
	for i = 1:length(od)
		od(i).LFP = [];
	end

	% specify interesting range
	uts = zeros(length(od(1).stimulus),1);
	uts(10e3:end-5e3) = true;
	for i = 1:length(od)
		od(i).use_this_segment = uts;
	end

	% back out all filters
	od = backOutFilters(od);
end

figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on
clear ax
ax(1) = subplot(2,2,1); hold on
ax(2) = subplot(2,2,2); hold on
for i = 3:6
	ax(i) = subplot(2,4,2+i); hold on
end

%  ######  ######## #### ##     ## ##     ## ##       ##     ##  ######  
% ##    ##    ##     ##  ###   ### ##     ## ##       ##     ## ##    ## 
% ##          ##     ##  #### #### ##     ## ##       ##     ## ##       
%  ######     ##     ##  ## ### ## ##     ## ##       ##     ##  ######  
%       ##    ##     ##  ##     ## ##     ## ##       ##     ##       ## 
% ##    ##    ##     ##  ##     ## ##     ## ##       ##     ## ##    ## 
%  ######     ##    #### ##     ##  #######  ########  #######   ######  


t = (1:length(od(1).stimulus))*1e-3;
plot(ax(1),t,nanmean([od.stimulus],2),'k')
xlabel(ax(1),'Time (s)')
ylabel(ax(1),'Stimulus (V)')
set(ax(1),'XLim',[0 70])

% ########  ########  ######  ########   #######  ##    ##  ######  ######## 
% ##     ## ##       ##    ## ##     ## ##     ## ###   ## ##    ## ##       
% ##     ## ##       ##       ##     ## ##     ## ####  ## ##       ##       
% ########  ######    ######  ########  ##     ## ## ## ##  ######  ######   
% ##   ##   ##             ## ##        ##     ## ##  ####       ## ##       
% ##    ##  ##       ##    ## ##        ##     ## ##   ### ##    ## ##       
% ##     ## ########  ######  ##         #######  ##    ##  ######  ######## 


for i = 1:length(od)
	plot(ax(2),t,nanmean([od(i).firing_rate],2),'Color',[0 0 0 0.1])
end
plot(ax(2),t,nanmean([od.firing_rate],2),'Color','k')
xlabel(ax(2),'Time (s)')
ylabel(ax(2),'Firing Rate (Hz)')
set(ax(2),'XLim',[0 70])


% ##     ##  ######         ##       #### ##    ## ########    ###    ########  
% ##     ## ##    ##        ##        ##  ###   ## ##         ## ##   ##     ## 
% ##     ## ##              ##        ##  ####  ## ##        ##   ##  ##     ## 
% ##     ##  ######         ##        ##  ## ## ## ######   ##     ## ########  
%  ##   ##        ##        ##        ##  ##  #### ##       ######### ##   ##   
%   ## ##   ##    ## ###    ##        ##  ##   ### ##       ##     ## ##    ##  
%    ###     ######  ###    ######## #### ##    ## ######## ##     ## ##     ## 


example_data = od(2);
% show response vs. linear projection; colour by mean stimulus in recent history window 
[~,excursions] = plotExcursions(example_data,ax(3),'data','firing_rate');

% make two time vectors, one defining when the stimulus is on, and one just for the whiffs
stim_on = false(length(example_data.stimulus),1);
stim_on(10e3:end-5e3) = true; 
whiff_times = false(length(example_data.stimulus),1);
for i = 1:length(excursions.ons)
	whiff_times(excursions.ons(i):excursions.offs(i)) = true;
end

% fit a NL just to the excursions
linear_model_data = ORNData;
linear_model_data.stimulus = nanmean(example_data.stimulus,2);
linear_model_data.firing_rate = nanmean(example_data.firing_rate,2);
linear_model_data.firing_projected = nanmean(example_data.firing_projected,2);
linear_model_data.use_this_segment = whiff_times;
[LN_model_data,ff] = fitNL(linear_model_data);

% show this best-fit NL on this plot
x = linear_model_data.firing_projected(linear_model_data.use_this_segment,:); x2 = x;
x2 = x2 - min(x2); x2 = x2/max(x2);
plot(ax(3),sort(x),max(linear_model_data.firing_rate(linear_model_data.use_this_segment,:))*ff(sort(x2)),'r')
xlabel(ax(3),'Proj. Stimulus (V)')
l = plot(ax(3),NaN,NaN);
linear_filter_r2 = rsquare(linear_model_data.firing_projected(whiff_times),linear_model_data.firing_rate(whiff_times));
legend(l,['r^2 = ' oval(linear_filter_r2)],'Location','southeast')

% ##     ##  ######         ##       ##    ##    ##     ##  #######  ########  ######## ##       
% ##     ## ##    ##        ##       ###   ##    ###   ### ##     ## ##     ## ##       ##       
% ##     ## ##              ##       ####  ##    #### #### ##     ## ##     ## ##       ##       
% ##     ##  ######         ##       ## ## ##    ## ### ## ##     ## ##     ## ######   ##       
%  ##   ##        ##        ##       ##  ####    ##     ## ##     ## ##     ## ##       ##       
%   ## ##   ##    ## ###    ##       ##   ###    ##     ## ##     ## ##     ## ##       ##       
%    ###     ######  ###    ######## ##    ##    ##     ##  #######  ########  ######## ######## 


% show response vs. LN model; colour by mean stimulus in recent history window 
uts = false(length(example_data.stimulus),1);
uts(10e3:end-5e3) = true;
LN_model_data.use_this_segment = uts;
plotExcursions(LN_model_data,ax(4),'data','firing_rate');
xlabel(ax(4),'LN Model Prediction (Hz)')
l = plot(ax(4),NaN,NaN);
LN_model_r2 = rsquare(LN_model_data.firing_projected(whiff_times),LN_model_data.firing_rate(whiff_times));
legend(l,['r^2 = ' oval(LN_model_r2)],'Location','southeast')

% ##     ##  ######     ########     ###       ##     ##  #######  ########  ######## ##       
% ##     ## ##    ##    ##     ##   ## ##      ###   ### ##     ## ##     ## ##       ##       
% ##     ## ##          ##     ##  ##   ##     #### #### ##     ## ##     ## ##       ##       
% ##     ##  ######     ##     ## ##     ##    ## ### ## ##     ## ##     ## ######   ##       
%  ##   ##        ##    ##     ## #########    ##     ## ##     ## ##     ## ##       ##       
%   ## ##   ##    ##    ##     ## ##     ##    ##     ## ##     ## ##     ## ##       ##       
%    ###     ######     ########  ##     ##    ##     ##  #######  ########  ######## ######## 


% fit a DA Model to the excursions
clear p
p.   s0 = -0.1164;
p.  n_z = 10.6250;
p.tau_z = 19.7499;
p.  n_y = 10.6250;
p.tau_y = 4.6377;
p.    C = 0.5848;
p.    A = 709.4439;
p.    B = 12.0094;
S = LN_model_data.stimulus;
[R,y,z,Ky,Kz] = DAModelv2(S,p);
DA_model_data = LN_model_data;
DA_model_data.firing_projected = R;
plotExcursions(DA_model_data,ax(5),'data','firing_rate');
l = plot(ax(5),NaN,NaN);
DA_model_r2 = rsquare(R(whiff_times),linear_model_data.firing_rate(whiff_times));
legend(l,['r^2 = ' oval(DA_model_r2)],'Location','southeast');
xlabel(ax(5),'DA Model Prediction (Hz)')


% now show that this is the key timescale of gain control
tau_gain = round(logspace(log10(50),4,50));
r2 = NaN*tau_gain;
for i = 1:length(tau_gain)
	p.tau_z = tau_gain(i)/p.n_z;
	R = DAModelv2(S,p);
	r2(i) = rsquare(R(whiff_times),linear_model_data.firing_rate(whiff_times));
end

% convert into fraction remaining variance explained
r2 = (r2 - linear_filter_r2)./(1-linear_filter_r2);


plot(ax(6),tau_gain,r2,'k+')
set(ax(6),'XScale','log','YLim',[0 1])
xlabel(ax(6),'Timescale of gain control (ms)')
ylabel(ax(6),['Remaining variance' char(10) 'explained by DA model'])

% show where the LN model is on this plot
plot(ax(6),tau_gain,(LN_model_r2-linear_filter_r2)/(1-linear_filter_r2)*(1+0*tau_gain),'r')


prettyFig('fs=18;','FixLogX=true;')

if being_published	
	snapnow	
	delete(gcf)
end


%% Controls
% In this section we establish controls for dynamical gain control. Negative controls are the LN model, which by definition cannot have gain control, and positive controls are the DA model, which does. 

figure('outerposition',[0 0 1500 800],'PaperUnits','points','PaperSize',[1500 800]); hold on
clear ax
for i = 8:-1:1
	ax(i) = subplot(2,4,i); hold on
end



%        ##     ## ########     ######   #######  ##    ## ######## ########   #######  ##       
%        ##     ## ##          ##    ## ##     ## ###   ##    ##    ##     ## ##     ## ##       
%        ##     ## ##          ##       ##     ## ####  ##    ##    ##     ## ##     ## ##       
% ###### ##     ## ######      ##       ##     ## ## ## ##    ##    ########  ##     ## ##       
%         ##   ##  ##          ##       ##     ## ##  ####    ##    ##   ##   ##     ## ##       
%          ## ##   ##          ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##       
%           ###    ########     ######   #######  ##    ##    ##    ##     ##  #######  ######## 


% use the LN model
control_data = ORNData;
q.   tau1 = 10;
q.   tau2 = 90.9375;
q.      A = 0.1484;
q.      n = 2;
q. offset = -0.0483;
q. Hill_A = 90.0547;
q.Hill_Kd = 0.6875;
q. Hill_n = 4;
[R,K] = pLNModel(linear_model_data.stimulus,q);
R = R + randn(length(R),1);
R(R<0) = 0;
control_data.stimulus = linear_model_data.stimulus;
control_data.firing_rate = R;
control_data.use_this_segment = stim_on;
control_data.regularisation_factor = 1e-1;
control_data.filtertime_firing = 1e-3*(1:length(K));
control_data.K_firing = K;


% show response vs. linear projection; colour by mean stimulus in recent history window 
[~,excursions] = plotExcursions(control_data,ax(1),'data','firing_rate');

% make two time vectors, one defining when the stimulus is on, and one just for the whiffs
stim_on = false(length(example_data.stimulus),1);
stim_on(10e3:end-5e3) = true; 
whiff_times = false(length(example_data.stimulus),1);
for i = 1:length(excursions.ons)
	whiff_times(excursions.ons(i):excursions.offs(i)) = true;
end

% fit a NL just to the excursions
control_data.use_this_segment = whiff_times;
[control_data_LN,ff] = fitNL(control_data);

% show this best-fit NL on this plot
x = control_data.firing_projected(control_data.use_this_segment,:); 
plot(ax(1),sort(x),max(control_data.firing_rate(control_data.use_this_segment,:))*ff(sort(x)),'r')
xlabel(ax(1),'Proj. Stimulus (V)')
l = plot(ax(1),NaN,NaN);
linear_filter_r2 = rsquare(control_data.firing_projected(control_data.use_this_segment),control_data.firing_rate(control_data.use_this_segment));
legend(l,['r^2 = ' oval(linear_filter_r2)],'Location','southeast')


% show response vs. LN model; colour by mean stimulus in recent history window 
uts = false(length(example_data.stimulus),1);
uts(10e3:end-5e3) = true;
control_data_LN.use_this_segment = uts;
plotExcursions(control_data_LN,ax(2),'data','firing_rate');
xlabel(ax(2),'LN Model Prediction (Hz)')
l = plot(ax(2),NaN,NaN);
LN_model_r2 = rsquare(control_data_LN.firing_projected(whiff_times),control_data_LN.firing_rate(whiff_times));
legend(l,['r^2 = ' oval(LN_model_r2)],'Location','southeast')

% fit a DA model
clear p
p. s0 = -0.1809;
p.  n_z = 2;
p.tau_z = 160.2500;
p.  n_y = 2;
p.tau_y = 8.8490;
p.    C = 1.3916e-06;
p.    A = 224.4745;
p.    B = 0.2726;
S = control_data_LN.stimulus;
R = DAModelv2(S,p);
DA_model_data = control_data_LN;
DA_model_data.firing_projected = R;
plotExcursions(DA_model_data,ax(3),'data','firing_rate');
l = plot(ax(3),NaN,NaN);
DA_model_r2 = rsquare(R(whiff_times),control_data_LN.firing_rate(whiff_times));
legend(l,['r^2 = ' oval(DA_model_r2)],'Location','southeast');
xlabel(ax(3),'DA Model Prediction (Hz)')


% now show that this is the key timescale of gain control
tau_gain = round(logspace(log10(50),4,50));
r2 = NaN*tau_gain;
for i = 1:length(tau_gain)
	p.tau_z = tau_gain(i)/p.n_z;
	R = DAModelv2(S,p);
	r2(i) = rsquare(R(whiff_times),control_data_LN.firing_rate(whiff_times));
end

% convert into fraction remaining variance explained
r2 = (r2 - linear_filter_r2)./(1-linear_filter_r2);


plot(ax(4),tau_gain,r2,'k+')
set(ax(4),'XScale','log','YLim',[0 1])
xlabel(ax(4),'Timescale of gain control (ms)')
ylabel(ax(4),['Remaining variance' char(10) 'explained by DA model'])

% show where the LN model is on this plot
plot(ax(4),tau_gain,(LN_model_r2-linear_filter_r2)/(1-linear_filter_r2)*(1+0*tau_gain),'r')



%        ##     ## ########     ######   #######  ##    ## ######## ########   #######  ##       
%   ##   ##     ## ##          ##    ## ##     ## ###   ##    ##    ##     ## ##     ## ##       
%   ##   ##     ## ##          ##       ##     ## ####  ##    ##    ##     ## ##     ## ##       
% ###### ##     ## ######      ##       ##     ## ## ## ##    ##    ########  ##     ## ##       
%   ##    ##   ##  ##          ##       ##     ## ##  ####    ##    ##   ##   ##     ## ##       
%   ##     ## ##   ##          ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##       
%           ###    ########     ######   #######  ##    ##    ##    ##     ##  #######  ######## 

% use the DA model
control_data = ORNData;
clear p
p.   s0 = -0.1164;
p.  n_z = 10.6250;
p.tau_z = 19.7499;
p.  n_y = 10.6250;
p.tau_y = 4.6377;
p.    C = 0.5848;
p.    A = 709.4439;
p.    B = 12.0094;
S = linear_model_data.stimulus;
R = DAModelv2(S,p);
R = R + randn(length(R),1);
R(R<0) = 0;
control_data.stimulus = linear_model_data.stimulus;
control_data.firing_rate = R;
control_data.use_this_segment = stim_on;
control_data.regularisation_factor = 1e-1;
control_data = backOutFilters(control_data);

% show response vs. linear projection; colour by mean stimulus in recent history window 
[~,excursions] = plotExcursions(control_data,ax(5),'data','firing_rate');

% make two time vectors, one defining when the stimulus is on, and one just for the whiffs
stim_on = false(length(example_data.stimulus),1);
stim_on(10e3:end-5e3) = true; 
whiff_times = false(length(example_data.stimulus),1);
for i = 1:length(excursions.ons)
	whiff_times(excursions.ons(i):excursions.offs(i)) = true;
end

% fit a NL just to the excursions
control_data.use_this_segment = whiff_times;
[control_data_LN,ff] = fitNL(control_data);

% show this best-fit NL on this plot
x = control_data.firing_projected(control_data.use_this_segment,:); 
plot(ax(5),sort(x),max(control_data.firing_rate(control_data.use_this_segment,:))*ff(sort(x)),'r')
xlabel(ax(5),'Proj. Stimulus (V)')
l = plot(ax(5),NaN,NaN);
linear_filter_r2 = rsquare(control_data.firing_projected(control_data.use_this_segment),control_data.firing_rate(control_data.use_this_segment));
legend(l,['r^2 = ' oval(linear_filter_r2)],'Location','southeast')


% show response vs. LN model; colour by mean stimulus in recent history window 
uts = false(length(example_data.stimulus),1);
uts(10e3:end-5e3) = true;
control_data_LN.use_this_segment = uts;
plotExcursions(control_data_LN,ax(6),'data','firing_rate');
xlabel(ax(6),'LN Model Prediction (Hz)')
l = plot(ax(6),NaN,NaN);
LN_model_r2 = rsquare(control_data_LN.firing_projected(whiff_times),control_data_LN.firing_rate(whiff_times));
legend(l,['r^2 = ' oval(LN_model_r2)],'Location','southeast')

% fit a DA model
clear p
p.   s0 = -0.1164;
p.  n_z = 2;
p.tau_z = 70.6249;
p.  n_y = 2;
p.tau_y = 20.8877;
p.    C = 0.4285;
p.    A = 694.5689;
p.    B = 12.0094;
S = control_data_LN.stimulus;
[R,y,z,Ky,Kz] = DAModelv2(S,p);
DA_model_data = control_data_LN;
DA_model_data.firing_projected = R;
plotExcursions(DA_model_data,ax(7),'data','firing_rate');
l = plot(ax(7),NaN,NaN);
DA_model_r2 = rsquare(R(whiff_times),control_data_LN.firing_rate(whiff_times));
legend(l,['r^2 = ' oval(DA_model_r2)],'Location','southeast');
xlabel(ax(7),'DA Model Prediction (Hz)')


% now show that this is the key timescale of gain control
tau_gain = round(logspace(log10(50),4,50));
r2 = NaN*tau_gain;
for i = 1:length(tau_gain)
	p.tau_z = tau_gain(i)/p.n_z;
	R = DAModelv2(S,p);
	r2(i) = rsquare(R(whiff_times),control_data_LN.firing_rate(whiff_times));
end

% convert into fraction remaining variance explained
r2 = (r2 - linear_filter_r2)./(1-linear_filter_r2);


plot(ax(8),tau_gain,r2,'k+')
set(ax(8),'XScale','log','YLim',[0 1])
xlabel(ax(8),'Timescale of gain control (ms)')
ylabel(ax(8),['Remaining variance' char(10) 'explained by DA model'])

% show where the LN model is on this plot
plot(ax(8),tau_gain,(LN_model_r2-linear_filter_r2)/(1-linear_filter_r2)*(1+0*tau_gain),'r')

for i = 1:4
	title(ax(i),'Synthetic Data: LN Model')
end
for i = 5:8
	title(ax(i),'Synthetic Data: DA Model')
end

prettyFig('fs=18;','FixLogX=true;')

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


