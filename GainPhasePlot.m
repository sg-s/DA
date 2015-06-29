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
	std_firing_rate = gain; 


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
		temp_std_firing_rate = [];
		for k = find(alldata(i).paradigm==j)
			% calculate stimulus variations

			y = alldata(i).resp(:,k);
			y = y(:);
			x = resp0(:); 
			ff = fit(x,y,'poly1');
			temp_gain = [temp_gain ff.p1];
			temp_mean_firing_rate = [temp_mean_firing_rate mean(y)];
			temp_std_firing_rate = [temp_std_firing_rate std(y)];


		end
		mean_firing_rate(j) = mean(temp_mean_firing_rate);


		stim_err (j) = std(temp_mean_firing_rate)/(sqrt(length(temp_mean_firing_rate)));

		gain(j) = mean(temp_gain);
		gain_err(j) = std(temp_gain)/(sqrt(length(temp_gain)));

		std_firing_rate(j) = mean(temp_std_firing_rate);
		



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
	end

	plot(gain,mean_firing_rate,'-+','Color',c)


end	

