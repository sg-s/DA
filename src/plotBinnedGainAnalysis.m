% plotBinnedGainAnalysis.m
% 
% created by Srinivas Gorur-Shandilya at 6:30 , 15 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [plot_handles,rho] = plotBinnedGainAnalysis(plot_here,s,p,r,plot_options)

plot_handles = [];

if plot_options.use_mean
	history_length = round(plot_options.history_length);
	K = ones(history_length,1);
	shat = filter(K,history_length,s);
else
	error('not coded')
end

% bin data
l = labelByPercentile(shat,plot_options.nbins);

stim_mean = NaN*(1:max(l));
stim_std = NaN*(1:max(l));
gain = NaN*(1:max(l));
gain_conf_intervals = NaN(2,max(l));
gain_fit_r2 = NaN*(1:max(l));

for i = 1:max(l)
	this_x = p(l==i);
	this_y = r(l==i);
	[ff,gof] = fit(this_x,this_y,'poly1');
	gain(i) = ff.p1;
	temp = confint(ff);
	gain_conf_intervals(:,i) = temp(:,1);
	gain_fit_r2(i) = gof.rsquare;
	stim_mean(i) = nanmean(shat(l==i));
	stim_std(i) = nanstd(shat(l==i));
end

keyboard

rm_this = gain<0;
gain(rm_this) = [];
stim_mean(rm_this) = [];
stim_std(rm_this) = [];
gain_conf_intervals(:,rm_this) = [];

if plot_options.make_plot
	plot_handles = errorbar(plot_here,stim_mean,gain,diff(gain_conf_intervals));
end

rho = spear(stim_mean(:),gain(:));