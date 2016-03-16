% fastGainControlAnalysis.m
% function that acts on a ORNData object and performs a fast gain analysis on it
% usage:
% fastGainControlAnalysis(ax,orn_data);
% where ax is a handle array pointing to 4 axes
% 
% you can get all the options using
% options = fastGainControlAnalysis;
% 
% created by Srinivas Gorur-Shandilya at 6:46 , 10 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [varargout] = fastGainControlAnalysis(ax,orn_data,varargin)

% get options from dependencies 
options = getOptionsFromDeps(mfilename);

% options and defaults
options.recompute_DA_fit = true;
options.max_exc_length = 500;
options.history_length = 300;
options.data = 'firing_rate'; 

if nargout && ~nargin 
	varargout{1} = options;
	return
end

% validate and accept options
if iseven(length(varargin))
	for ii = 1:2:length(varargin)-1
	temp = varargin{ii};
    if ischar(temp)
    	if ~any(find(strcmp(temp,fieldnames(options))))
    		disp(['Unknown option: ' temp])
    		disp('The allowed options are:')
    		disp(fieldnames(options))
    		error('UNKNOWN OPTION')
    	else
    		options = setfield(options,temp,varargin{ii+1});
    	end
    end
end
elseif isstruct(varargin{1})
	% should be OK...
	options = varargin{1};
else
	error('Inputs need to be name value pairs')
end

assert(orn_data.n_trials == 1,'fastGainControlAnalysis expects only one trial of data, which is averaged over all trials')

% grab the core region where the stim is on
stim_on = orn_data.use_this_segment;



if strcmp(options.data,'LFP')
	% make sure the LFP is uniformly negative
	orn_data.LFP = orn_data.LFP - max(orn_data.LFP(stim_on));
end

orn_data = backOutFilters(orn_data);

% show response vs. linear projection; colour by mean stimulus in recent history window 
[~,excursions] = plotExcursions(orn_data,ax(1),'max_exc_length',options.max_exc_length,'history_length',options.history_length,'data',options.data);

assert(length(excursions.ons)>25,'Too few excursions!')

% make two time vectors, one defining when the stimulus is on, and one just for the whiffs
whiff_times = false(length(orn_data.stimulus),1);
for i = 1:length(excursions.ons)
	whiff_times(excursions.ons(i):excursions.offs(i)) = true;
end

% fit a NL just to the excursions
orn_data.use_this_segment = whiff_times;
[orn_data_LN,ff] = fitNL(orn_data,options.data);

% show this best-fit NL on this plot
if strcmp(options.data,'LFP')
	x = -orn_data.LFP_projected(orn_data.use_this_segment,:); 
	plot(ax(1),sort(x), -min(orn_data.LFP(orn_data.use_this_segment,:))*ff(sort(x)),'r')
else
	x = orn_data.firing_projected(orn_data.use_this_segment,:); 
	plot(ax(1),sort(x),max(orn_data.firing_rate(orn_data.use_this_segment,:))*ff(sort(x)),'r')
end


l = plot(ax(1),NaN,NaN);
if strcmp(options.data,'LFP')
	xlabel(ax(1),'-Proj. Stimulus (V)')
	linear_filter_r2 = rsquare(orn_data.LFP_projected(orn_data.use_this_segment),orn_data.LFP(orn_data.use_this_segment));
else
	xlabel(ax(1),'Proj. Stimulus (V)')
	linear_filter_r2 = rsquare(orn_data.firing_projected(orn_data.use_this_segment),orn_data.firing_rate(orn_data.use_this_segment));
end
legend(l,['r^2 = ' oval(linear_filter_r2)],'Location','southeast')


% show response vs. LN model; colour by mean stimulus in recent history window 
orn_data_LN.use_this_segment = stim_on;
plotExcursions(orn_data_LN,ax(2),'data',options.data,'max_exc_length',options.max_exc_length);
l = plot(ax(2),NaN,NaN);
if strcmp(options.data,'LFP')
	xlabel(ax(2),'-LN Model Prediction (mV)')
	LN_model_r2 = rsquare(orn_data_LN.LFP_projected(whiff_times),orn_data_LN.LFP(whiff_times));
else
	xlabel(ax(2),'LN Model Prediction (Hz)');
	LN_model_r2 = rsquare(orn_data_LN.firing_projected(whiff_times),orn_data_LN.firing_rate(whiff_times));
end

legend(l,['r^2 = ' oval(LN_model_r2)],'Location','southeast')

% fit a DA model
clear d p0
p0.   s0 = -0.1164;
p0.  n_z = 10.6250;
p0.tau_z = 19.7499;
p0.  n_y = 10.6250;
p0.tau_y = 4.6377;
p0.    C = 0.5848;
p0.    A = 709.4439;
p0.    B = 12.0094;

d.stimulus = orn_data.stimulus;
if strcmp(options.data,'LFP')
	d.response = -orn_data.LFP;
else
	d.response = orn_data.firing_rate;
end
d.response(~orn_data.use_this_segment) = NaN;
p = cache(dataHash(d));

if isempty(p)
	p = fitModel2Data(@DAModelv2,d,'p0',p0);
	cache(dataHash(d),p);
else
	if options.recompute_DA_fit
		cache(dataHash(d),[]);
		p = fitModel2Data(@DAModelv2,d,'p0',p);
		cache(dataHash(d),p);
	end
end

S = orn_data_LN.stimulus;
R = DAModelv2(S,p);
DA_model_data = orn_data_LN;
if strcmp(options.data,'LFP')
	DA_model_data.LFP_projected = -R;
else
	DA_model_data.firing_projected = R;
end
plotExcursions(DA_model_data,ax(3),'data',options.data,'max_exc_length',options.max_exc_length);
l = plot(ax(3),NaN,NaN);
if strcmp(options.data,'LFP')
	xlabel(ax(3),'-DA Model Prediction (mV)')
	DA_model_r2 = rsquare(DA_model_data.LFP_projected(whiff_times),DA_model_data.LFP(whiff_times));
else
	xlabel(ax(3),'DA Model Prediction (Hz)');
	DA_model_r2 = rsquare(DA_model_data.firing_projected(whiff_times),DA_model_data.firing_rate(whiff_times));
end
legend(l,['r^2 = ' oval(DA_model_r2)],'Location','southeast');


% now show that this is the key timescale of gain control
tau_gain = round(logspace(log10(50),4,50));
r2 = NaN*tau_gain;
for i = 1:length(tau_gain)
	p.tau_z = tau_gain(i)/p.n_z;
	R = DAModelv2(S,p);
	r2(i) = rsquare(R(whiff_times),orn_data_LN.firing_rate(whiff_times));
end

% convert into fraction remaining variance explained
r2 = (r2 - linear_filter_r2)./(1-linear_filter_r2);
[~,loc] = max(r2);


r2_plot_data.l = plot(ax(4),tau_gain,r2,'k+');
set(ax(4),'XScale','log','YLim',[0 1])
xlabel(ax(4),'Timescale of gain control (ms)')
ylabel(ax(4),['Remaining variance' char(10) 'explained by DA model'])

% show where the LN model is on this plot
plot(ax(4),tau_gain,(LN_model_r2-linear_filter_r2)/(1-linear_filter_r2)*(1+0*tau_gain),'r')
r2_plot_data.peak_tau = tau_gain(loc);

varargout{1} = r2_plot_data;