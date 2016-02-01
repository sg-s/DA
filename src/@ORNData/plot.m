%ORNDATA/plot.m
% 
% created by Srinivas Gorur-Shandilya at 9:20 , 01 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [plot_handles] = plot(varargin)

% options and defaults
showr2 = false;
ioCurve_type = 'pwlinear'; % can also be "dots"
nbins = 30;
plot_type = 'trial-wise'; % can also be "trial-wise", or "mean"
grouping = [];
plot_here = [];

% defensive programming
assert(length(varargin)>1,'Not enough input arguments.')
o = varargin{1};
assert(isa(o,'ORNData'),'1st argument should be an object from the ORNData class')
varargin(1) = [];

% figure out where to plot
temp = varargin{1};
if strcmp(class(temp),'matlab.graphics.axis.Axes')
	% user has supplied a handlle to a axis; make sure it is valid
	if min(isvalid(temp))
		plot_here = varargin{1};

	else
	end
	varargin(1) = [];
else
end

% figure out WHAT to plot
allowed_plots = {'L.','muSigma.','ioCurve.','LN.','weber.','laughlin.','gainAnalysis.'};
temp = varargin{1};
if strcmp(class(temp),'char')
	% user has supplied a string, make sure we can understand it 
	ok = false;
	for i = 1:length(allowed_plots)
		if ~isempty(strfind(varargin{1},allowed_plots{i}))
			ok = true;
		end
	end
	if ok 
		plot_what = varargin{1};
		varargin(1) = [];
	else
		disp('I dont know what you want me to plot. Allowed plot arguments are:')
		disp(allowed_plots')
		error('Invalid plot argument')
	end
else
	error('Expected a string specifying what to plot, got something else instead.')
end

% handle any additional arguments (they have to be options in name value syntax )
if length(varargin) 
	if iseven(length(varargin))
		for ii = 1:2:length(varargin)-1
			temp = varargin{ii};
        	if ischar(temp)
        		eval(strcat(temp,'=varargin{ii+1};'));
        	end
    	end
	else
    	error('Inputs need to be name value pairs')
    end
end

% respect use_this_segment and *se_these_trials properties
if isempty(o.use_these_trials)
	utt = 1:o.n_trials;
else
	utt = o.use_these_trials;
end
if isempty(o.use_this_segment)
	uts = 1:length(o.stimulus);
else
	uts = o.use_this_segment;
end


if strfind(plot_what,'ioCurve.') 
	if isempty(plot_here)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end
	% only plot nonlinearity. figure out what to plot
	if strfind(plot_what,'firing')
		pred = o.firing_projected(uts,utt);
		resp = o.firing_rate(uts,utt);

	elseif strfind(plot_what,'LFP')
		error('73 not coded')
	else
		error('What to plot not specified. You told me to plot the ioCurve, but the only things I can plot the ioCurve of are "firing" or "LFP"')
	end

	clear plot_options
	plot_options.ioCurve_type = ioCurve_type;
	plot_options.nbins = nbins;
	plot_options.showr2 = showr2;
	plot_options.plot_type = plot_type;
	plot_handles = plotNonlinearity(plot_here,pred,resp,grouping,plot_options);
end
