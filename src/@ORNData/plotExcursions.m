% @ORNData/plotExcursions.m
% plots excursions (defined as >10% peak responses) for a ORNData class object
% usage
% plotExcursions(od) % makes a plot showing firing rate excursions, colours excursions
% plotExcursions(plot_here,od) % where plot_here is a handle to a axes
% plotExcursions(plot_here,od,'history_length',500) % in ms
% plotExcursions(plot_here,od,'data','LFP')
% plotExcursions(plot_here,od,'data','firing_rate')
% plotExcursions(plot_here,od,'normalise',true)
% created by Srinivas Gorur-Shandilya at 4:04 , 07 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [varargout] = plotExcursions(varargin)

% options and defaults
options.history_length = 500;
options.excursion_thresh = .1;
options.min_exc_length = 50;
options.max_exc_length = 500;
options.data = 'firing_rate';
options.min_excursion_r2 = .8;

if nargout && ~nargin
	varargout{1} = options;
end

% grab the ORNData object from the arguments
for i = 1:length(varargin)
	if isa((varargin{i}),'ORNData')
		o = varargin{i};
		varargin(i) = [];
		break
	end
end

% grab the axes object, if it exists, from the arguments
plot_here = [];
for i = 1:length(varargin)
	if isa((varargin{i}),'matlab.graphics.axis.Axes')
		plot_here = varargin{i};
		varargin(i) = [];
		break
	end
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
else
	error('Inputs need to be name value pairs')
end

% make sure axes handle is valid
if ~isempty(plot_here)
	if ~isvalid(plot_here)
		plot_here = [];
	end
end

if isempty(plot_here)
	clear plot_here % so that plot_here gets the right class (axes object)
	figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
	plot_here = gca;
	hold(plot_here,'on')
end


% respect use_this_segment and use_these_trials properties
if isempty(o.use_these_trials)
	utt = true(o.n_trials,1);
else
	utt = o.use_these_trials;
	if ~isempty(grouping)
		grouping = grouping(utt);
	end
	
end
if isempty(o.use_this_segment)
	uts = true(length(o.stimulus),1);
else
	uts = o.use_this_segment;
end


stim = nanmean(o.stimulus(uts,utt),2);

% excursions are always defined in the firing rate space. 
f = nanmean(o.firing_rate(uts,utt),2);

% find all the excursions in the firing rate
f = f-min(f);
f = f/max(f);
[ons,offs] = computeOnsOffs(f>options.excursion_thresh);

% excursions should be a certain length
rm_this = (offs-ons) < options.min_exc_length | (offs-ons) > options.max_exc_length;
ons(rm_this) = []; offs(rm_this) = [];

% figure out what to plot
if strfind(options.data,'firing_rate')
	pred = nanmean(o.firing_projected(uts,utt),2);
	resp = nanmean(o.firing_rate(uts,utt),2);
	ylabel(plot_here(1),'Firing rate (Hz)')

	% throw out some junk points
	rm_this = pred > 2*max(resp);
	pred(rm_this) = []; resp(rm_this) = [];

elseif strfind(options.data,'LFP')
	pred = -nanmean(o.LFP_projected(uts,utt),2);
	resp = -nanmean(o.LFP(uts,utt),2);
	ylabel(plot_here(1),'-LFP (mV)')
else
	error('I dont understand what data to plot.')
end

% find the gains in all windows and also grab the data to plot
[gain,gain_err,plot_data] = findGainInWindows(ons,offs,pred,resp);

% throw out some crappy data
rm_this = gain < 0 | gain_err < options.min_excursion_r2;
gain(rm_this) = [];
gain_err(rm_this) = [];
plot_data(rm_this) = [];
ons(rm_this) = [];
offs(rm_this) = [];


% show the example history length
temp = computeSmoothedStimulus(stim,round(options.history_length));
shat = NaN*ons;
for j = 1:length(ons)
	shat(j) = mean(temp(ons(j):offs(j)));
end


% colour the plot 
shat = shat - nanmin(shat);
shat = shat/nanmax(shat);
c = parula(100);
for i = length(ons):-1:1
	try
		plot_handles(i) = plot(plot_here,pred(ons(i):offs(i)),resp(ons(i):offs(i)),'.','Color',c(1+floor(shat(i)*99),:));
	catch
	end
end
	

% also return the excursions.
excursions.ons = ons + find(uts,1,'first');
excursions.offs = offs + find(uts,1,'first');

varargout{1} = plot_handles;
varargout{2} = excursions;
