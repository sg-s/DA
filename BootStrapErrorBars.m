% BootStrap Slopes.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
% x is a structure with the following fields:
% x.data -- the data
% x.prediction -- the prediction
% x.stimulus -- the stimulus
% x.time -- a time vector (uniformly spaced)
% all these vectors are equally long
function [low_slopes, high_slopes, p] = BootStrapErrorBars(x,history_lengths,fraction)

% unpack data
f = x.response(:);
fp = x.prediction(:);
stimulus = x.stimulus(:);
t = x.time(:);

% figure out sampling rate
dt = mean(diff(t));
hl = round(history_lengths/dt);


% compute shat
shat = ComputeSmoothedStimulus(stimulus,hl);

% parameters
nrep = 100; % how many times do we bootstrap the data?

% initialise outputs
p = NaN(1,length(hl));  % stores p-values for each history length

low_slopes.bootstrap = NaN(nrep,length(hl));
low_slopes.data = NaN(1,length(hl));

high_slopes.bootstrap= NaN(nrep,length(hl));
high_slopes.data = NaN(1,length(hl));

n = floor(fraction*length(shat));

for i = 1:length(hl) % for each history length
	textbar(i,length(hl))
	% do low slopes
	this_shat = shat(i,:);
	this_shat(1:hl(i)) = Inf; % the initial segment where we can't estimate shat is excluded
	this_shat(isnan(this_shat)) = Inf; % exclude NaNs
	[~, idx] = sort(this_shat,'ascend');
	idx = idx(:);

	% calculate the slopes
	f_low = f(idx(1:n));
	fp_low = fp(idx(1:n));

	% strip NaN
	f_low(isnan(fp_low)) = [];
	fp_low(isnan(fp_low)) = [];

	% fit lines
	flow = fit(fp_low,f_low,'Poly1');
	low_slopes.data(i) = flow.p1;


	shift =  randi(round(length(idx)),nrep,1);
	for j = 1:nrep  % bootstrap this nrep times
		shifted_idx = circshift(idx,[shift(j) 0]);
		f_low = f(shifted_idx(1:n));
		fp_low = fp(shifted_idx(1:n));

		% strip NaN
		f_low(isnan(fp_low)) = [];
		fp_low(isnan(fp_low)) = [];

		% fit lines
		flow = fit(fp_low,f_low,'Poly1');
		low_slopes.bootstrap(j,i) = flow.p1;
	end
	clear j





	% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	% do the same thing for the high slopes
	this_shat = shat(i,:);
	this_shat(1:hl(i)) = -Inf; % the initial segment where we can't estimate shat is excluded
	this_shat(isnan(this_shat)) = -Inf; % exclude NaNs
	[~, idx] = sort(this_shat,'descend');
	idx = idx(:);

	% calculate the slopes
	f_low = f(idx(1:n));
	fp_low = fp(idx(1:n));

	% strip NaN
	f_low(isnan(fp_low)) = [];
	fp_low(isnan(fp_low)) = [];

	% fit lines
	flow = fit(fp_low,f_low,'Poly1');
	high_slopes.data(i) = flow.p1;

	for j = 1:nrep  % bootstrap this many times
		shift =  randi(length(idx));
		shifted_idx = circshift(idx,[shift 0]);
		f_low = f(shifted_idx(1:n));
		fp_low = fp(shifted_idx(1:n));

		% strip NaN
		f_low(isnan(fp_low)) = [];
		fp_low(isnan(fp_low)) = [];

		% fit lines
		flow = fit(fp_low,f_low,'Poly1');
		high_slopes.bootstrap(j,i) = flow.p1;
	end
	clear j


	% find the absolute value of the difference between high and low slopes
	a = abs(low_slopes.bootstrap(:,i) - high_slopes.bootstrap(:,i));
	a0 = abs(low_slopes.data(i) - high_slopes.data(i));

	% calculate p 
	p(i) = sum(a>a0)/nrep;

end
clear i

