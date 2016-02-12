% plotInstGainVsStim.m
% 
% created by Srinivas Gorur-Shandilya at 3:32 , 02 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [plot_handles,rho] = plotInstGainVsStim(plot_here,inst_gain,inst_gain_err,stim,resp,plot_options)

if plot_options.use_mean
	history_length = round(plot_options.history_length);
	K = ones(history_length,1);
	stim = filter(K,history_length,stim);
else
	error('not coded')
end

% throw out some data
if min(resp) >= 0
	% it's firing rate
	rm_this = (inst_gain_err<plot_options.min_inst_gain_r2 | inst_gain < 0 | isinf(inst_gain) | resp < plot_options.min_inst_gain_firing);
else
	% it's LFP
	rm_this = (inst_gain_err<plot_options.min_inst_gain_r2 | inst_gain < 0 | isinf(inst_gain));

end


rm_this(1:history_length*2) = true;
inst_gain(rm_this) = [];
stim(rm_this) = [];


if plot_options.make_plot
	if strcmp(plot_options.data_bin_type,'pwlinear')
		axes(plot_here);
		plot_handles = plotPieceWiseLinear(stim,inst_gain,'nbins',plot_options.nbins,'Color',plot_options.colour,'use_std',plot_options.use_std);
	elseif strcmp(plot_options.data_bin_type,'dots')
		plot_handles = plot(plot_here,stim(1:plot_options.nbins:end),inst_gain(1:plot_options.nbins:end),'Color',plot_options.colour,'Marker','.','LineStyle','none');
	end
		
else
	plot_handles = [];
end

if length(inst_gain) > 1e3
	rho = spear(stim(1:10:end),inst_gain(1:10:end));
else
	rho = spear(stim,inst_gain);
end

% also fit a power law and get the exponent
% ff = fit(stim,inst_gain,'power1');
% weber_exponent = ff.b;

