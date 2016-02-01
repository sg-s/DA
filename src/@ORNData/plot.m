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
history_length = 500; % in ms, or the time step of the data (which is the same)
hl_min = 200; % in ms
hl_max = 1e4;
history_lengths = unique([logspace(log10(hl_min),log10(500),15) logspace(log10(500),log10(hl_max),15)]);


% defensive programming
assert(length(varargin)>1,'Not enough input arguments.')
o = varargin{1};
assert(isa(o,'ORNData'),'1st argument should be an object from the ORNData class')
varargin(1) = [];

% figure out WHERE to plot
temp = varargin{1};
if isa(temp,'matlab.graphics.axis.Axes')
	% user has supplied a handlle to a axis; make sure it is valid
	if min(isvalid(temp))
		plot_here = varargin{1};

	else
	end
	varargin(1) = [];
else
end

% figure out WHAT to plot
allowed_plots = {'pdf','Filter.','muSigma.','ioCurve.','LN.','weber.','laughlin.','gainAnalysis.','muSigmaR2.'};
temp = varargin{1};
if isa(temp,'char')
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
if ~isempty(varargin) 
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


if strfind(plot_what,'Filter.')
	if isempty(plot_here)
		clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end
	% only plot fitlers. figure out which filters to plot
	if strfind(plot_what,'firing')
		filtertime = o.filtertime_firing;
		K = o.K_firing;
	elseif strfind(plot_what,'LFP')
		filtertime = o.filtertime_LFP;
		K = o.K_LFP;
	else
		error('What to plot not specified. You told me to plot the ioCurve, but the only things I can plot the ioCurve of are "firing" or "LFP"')
	end

	if strcmp(plot_type,'trial-wise')
		plot_handles = plot(plot_here,filtertime,K,'Color',[.5 .5 .5]);
		plot_handles(end+1) = plot(plot_here,filtertime,nanmean(K,2),'Color','k','LineWidth',3);
	elseif strcmp(plot_type,'sem')
		axis(plot_here);
		plot_handles = errorShade(filtertime,nanmean(K,2),sem(K'),'Color',[0 0 0]);
	elseif strcmp(plot_type,'mean')
		plot_handles = plot(plot_here,filtertime,nanmean(K,2),'Color','k','LineWidth',2);
	else
		error('Unknown plot type')
	end
	xlabel(plot_here,'Filter Lag (s)')
	ylabel(plot_here,'Filter')

end

if strfind(plot_what,'ioCurve.') 
	if isempty(plot_here)
		clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end
	% only plot nonlinearity. figure out what to plot
	if strfind(plot_what,'firing')
		pred = o.firing_projected(uts,utt);
		resp = o.firing_rate(uts,utt);
		ylabel(plot_here,'Firing Rate (Hz)')

	elseif strfind(plot_what,'LFP')
			ylabel(plot_here,'\DeltaLFP (mV)')
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
	xlabel(plot_here,'Projected Stimulus (V)')

end

if strfind(plot_what,'LN.')
	 if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on 			
		plot_here(1) = subplot(1,2,1); hold on
		plot_here(2) = subplot(1,2,2); hold on
	end
	plot(o,plot_here(1),strrep(plot_what,'LN.','Filter.'));
	plot(o,plot_here(2),strrep(plot_what,'LN.','ioCurve.'));

end

if strfind(plot_what,'pdf.')
	% show a distribution of ... something
	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end

	if strfind(plot_what,'stimulus')
		xlabel('Stimulus (V)')
		x = o.stimulus(uts,utt);
	elseif strfind(plot_what,'firing') && ~strfind(plot_what,'firing_projected')
		xlabel('Firing Rate (Hz)')
		x = o.firing_rate(uts,utt);
	elseif strfind(plot_what,'firing_projected')
		xlabel('Projected Stimulus (V)')
		x = o.firing_projected(uts,utt);
	elseif strfind(plot_what,'LFP') && ~strfind(plot_what,'LFP_projected')
		xlabel('\DeltaLFP (mV)')
		x = o.LFP(uts,utt);
	elseif strfind(plot_what,'LFP_projected')
		x = o.LFP_projected(uts,utt);
		xlabel('Projected Stimulus (V)')
	elseif strfind(plot_what,'inst_gain')
		error('185 not coded')
	elseif strfind(plot_what,'inst_gain_r2')
		error('187 not coded')
	end

	ylabel(plot_here,'p.d.f')

	clear plot_options
	plot_options.ioCurve_type = ioCurve_type;
	plot_options.nbins = nbins;
	plot_options.showr2 = showr2;
	plot_options.plot_type = plot_type;
	plot_handles = plotPDF(plot_here,x,grouping,plot_options);
end

if any(strfind(plot_what,'muSigma.')) && ~any(strfind(plot_what,'muSigmaR2.'))
	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end

	if strfind(plot_what,'stimulus') 
		plot_this = o.stimulus(uts,utt);
		xlabel(plot_here,'\mu_{Stimulus} (V)')
		ylabel(plot_here,'\sigma_{Stimulus} (V)')
	elseif any(strfind(plot_what,'firing')) && ~any(strfind(plot_what,'firing_projected'))
		plot_this = o.firing_rate(uts,utt);
		xlabel(plot_here,'\mu_{Firing Rate} (Hz)')
		ylabel(plot_here,'\sigma_{Firing Rate} (Hz)')
	elseif strfind(plot_what,'firing_projected')
		plot_this = o.firing_projected(uts,utt);
		xlabel(plot_here,'\mu_{Projected Stim.} (V)')
		ylabel(plot_here,'\sigma_{Projected Stim.} (V)')
	elseif strfind(plot_what,'LFP') && ~strfind(plot_what,'LFP_projected')
		plot_this = o.LFP(uts,utt);
	elseif strfind(plot_what,'LFP_projected')
		plot_this = o.LFP_projected(uts,utt);
	end
	plot_this = plot_this(:);

	plot_handles = plotMuSigma(plot_here,plot_this,history_length,true);
end

if any(strfind(plot_what,'muSigmaR2.'))
	if isempty(plot_here)
	 	clear plot_here % so that plot_here gets the right class (axes object)
		figure('outerposition',[0 0 600 600],'PaperUnits','points','PaperSize',[1000 600]); hold on 			
		plot_here = gca;
	end

	if strfind(plot_what,'stimulus') 
		plot_this = o.stimulus(uts,utt);
	elseif any(strfind(plot_what,'firing')) && ~any(strfind(plot_what,'firing_projected'))
		plot_this = o.firing_rate(uts,utt);
	elseif strfind(plot_what,'firing_projected')
		plot_this = o.firing_projected(uts,utt);
	elseif strfind(plot_what,'LFP') && ~strfind(plot_what,'LFP_projected')
		plot_this = o.LFP(uts,utt);
	elseif strfind(plot_what,'LFP_projected')
		plot_this = o.LFP_projected(uts,utt);
	end
	ylabel(plot_here,'\rho')
	xlabel(plot_here,'History Length (ms)')
	plot_this = plot_this(:);

	r2 = NaN*history_lengths;
	for i = 1:length(history_lengths)
		[~,r2(i)] = plotMuSigma([],plot_this,history_lengths(i),false);
	end
	plot(plot_here,history_lengths,r2,'k+')
	set(plot_here,'XScale','log','YLim',[-1 1])

end






