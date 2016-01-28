% plotSpearmanRhoVsHistoryLengths.m
% this subfunction is called by the overloaded plot method in the ORNData class
% 
% created by Srinivas Gorur-Shandilya at 10:44 , 28 January 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function plot_handles = plotSpearmanRhoVsHistoryLengths(plot_here,o,mean_stim_bins,normalise_gain,normalise_preceding_stim,history_lengths,use_LFP)

if ~exist('use_LFP','var')
	use_LFP = false;
end	


stimulus = nanmean(o.stimulus,2);
if use_LFP
	prediction = nanmean(o.LFP_projected,2);
	response = nanmean(o.LFP,2);
else
	prediction = nanmean(o.firing_projected,2);
	response = nanmean(o.firing_rate,2);
end

% check the cache
temp.stimulus = stimulus;
temp.mean_stim_bins = mean_stim_bins;
temp.normalise_gain = normalise_gain;
temp.normalise_preceding_stim = normalise_preceding_stim;
temp.history_lengths = history_lengths;
temp.dt = o.dt;
hash = dataHash(temp);

rho = cache(hash);


if isempty(rho)
	% ok, calculate

	% throw out the first 5 seconds
	stimulus = stimulus(5e3:end);
	response = response(5e3:end);
	prediction = prediction(5e3:end);

	if normalise_gain
		% fix the gain to be exactly 1
		x = prediction(:);
		y = response(:);
		rm_this = isnan(x) | isnan(y);
		x(rm_this) = [];
		y(rm_this) = [];
		temp = fit(x,y,'poly1');
		prediction = prediction*temp.p1;
	end

	rho = NaN*history_lengths;


	for j = 1:length(history_lengths)
		hl = floor(history_lengths(j)/o.dt);
		shat = filter(ones(hl,1),hl,stimulus);

		shat(1:hl) = [];  % we can't estimate shat here 
		this_prediction = prediction; this_prediction(1:hl) = [];
		this_response = response; this_response(1:hl) = [];

		if normalise_preceding_stim
			shat = shat - min(shat);
			shat = shat/nanmax(shat);
		end


		l = labelByPercentile(shat,mean_stim_bins);
		x = NaN(mean_stim_bins,1);
		y = NaN(mean_stim_bins,1);

		for bi = 1:mean_stim_bins
			x(bi) = nanmean(shat(l==bi));
			p = (this_prediction(l==bi));
			r = (this_response(l==bi));
			rm_this = isnan(p) | isnan(r);
			p(rm_this) = []; r(rm_this) = [];
			[temp] = pca([p r]);
			y(bi) = temp(2,1)/temp(1,1);
			
		end
		rho(j) = spear(x,y);


	end
	% cache for later use
	cache(hash,rho)
end

plot_handles = plot(plot_here,history_lengths,rho,'k+');


