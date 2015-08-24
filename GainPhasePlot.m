% makes a phase plot of mean firing rate vs. slope for experiments where we use light and odor, one flickering and one constant. 
% 
% created by Srinivas Gorur-Shandilya at 10:34 , 01 May 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [max_f, max_gain] = GainPhasePlot(alldata,c,norm_f)


max_f = NaN(length(alldata),1);
max_gain = NaN(length(alldata),1);

for i = 1:length(alldata)

	% find which part of the trace we can use
	rm_this=(any(isnan(alldata(i).resp')));
	alldata(i).resp(rm_this,:) = [];
	alldata(i).stim(rm_this,:) = [];

	gain = NaN(length(unique(alldata(i).paradigm)),1);
	gain_err = gain;
	mean_firing_rate = gain;
	mean_firing_rate_err = gain; 


	% first grab some data about the no background case (so either odor flicker alone or light flicker alone)
	resp0= alldata(i).resp(:,alldata(i).paradigm==1);
	stim0= alldata(i).stim(:,alldata(i).paradigm==1);
	if width(resp0) > 1
		resp0 = mean2(resp0);
		stim0 = mean2(stim0);
	end

	
	for j = 1:length(alldata(i).paradigm) % for each paradigm
		temp_gain = [];
		temp_mean_firing_rate = [];
		ci = [];
		for k = find(alldata(i).paradigm==j)
			y = alldata(i).resp(:,k);
			y = y(:);
			x = resp0(:); 
			ff = fit(x,y,'poly1');
			temp_gain = [temp_gain ff.p1];
			ci_temp = confint(ff);
			ci = [ci ci_temp(1,1)];
			temp_mean_firing_rate = [temp_mean_firing_rate mean(y)];


		end

		if length(temp_mean_firing_rate) == 1
			gain(j) = temp_gain;
			gain_err(j) = temp_gain - ci;
			
			mean_firing_rate(j) = temp_mean_firing_rate;
			mean_firing_rate_err(j) = std(y);

		else

			mean_firing_rate(j) = mean(temp_mean_firing_rate);
			mean_firing_rate_err(j) = sem(temp_mean_firing_rate);

			gain(j) = sum((1./ci.*temp_gain)/(sum(1./ci))); % weighted by error in each fit
			gain_err(j) = sem(temp_gain);


		end
	end

	% save for output
	[max_f(i),idx] = max(mean_firing_rate);
	max_gain(i) = gain(idx);

	xt = {};
	for j = 1:length(alldata(i).ParadigmNames)
		t = alldata(i).ParadigmNames{j};
		a = strfind(t,'+');
		z = strfind(t,'V');
		xt{j} = t(a+1:z-1);
	end

	% normalise firing rates
	if norm_f
		mean_firing_rate = mean_firing_rate/mean_firing_rate(1);
		mean_firing_rate_err = mean_firing_rate_err/mean_firing_rate(1);
	end

	% also plot error bars
	errorbarxy(gain,mean_firing_rate,gain_err,mean_firing_rate_err,c)


end	

