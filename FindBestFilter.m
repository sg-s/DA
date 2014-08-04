% FindBestFilter.m
% created by Srinivas Gorur-Shandilya at 17:31 , 15 January 2014. Contact me at http://srinivas.gs/contact/
% finds the best filter, searching through the space of the free parameter in my and damon's filter estimation functions
function [K, diagnostics, filtertime] = FindBestFilter(stim,response,OnlyThesePoints,varargin)

% ensure inputs OK
stim = stim(:);
response = response(:);

% defaults
regmax = 100;
regmin = 1e-5;
algo = 1;
reg = [];
regtype=2;
filter_length = 333;
min_cutoff = -Inf;
offset = floor(filter_length/10);

if nargin < 3
	OnlyThesePoints = [];
end



% evaluate optional inputs
for i = 1:nargin-3
	eval(varargin{i})
end


% offset the stimulus and response a little bit to account for acausal filters (I know, weird)
stim = stim(offset:end);
response = response(1:end-offset+1);
filtertime = [-offset+1:filter_length-offset+1];

if algo == 1
	%% Chichilnisky's method. 


	if regmax == regmin
		% just use what's given, don't optimise
		flstr = strcat('filter_length=',mat2str(filter_length),';');
		regstr = strcat('reg=',mat2str(regmax),';');
		regstr2 = strcat('regtype=',mat2str(regtype),';');
		K = FitFilter2Data(stim,response,OnlyThesePoints,flstr,regstr,regstr2);
		diagnostics = [];
		return
	end

	nsteps = 20;
	ss = (log(regmax)-log(regmin))/nsteps;
	reg = log(regmin):ss:log(regmax);
	reg = exp(reg);
	if length(reg) > nsteps
		reg(end) = [];
	end
	 
	err = NaN(1,nsteps);
	filter_height = NaN(1,nsteps);
	filter_sum = NaN(1,nsteps);
	slope = NaN(1,nsteps);
	mcond = NaN(1,nsteps);


	for i = 1:nsteps
		% find the filter
		flstr = strcat('filter_length=',mat2str(filter_length),';');
		regstr = strcat('reg=',mat2str(reg(i)),';');
		regstr2 = strcat('regtype=',mat2str(regtype),';');
		[K, C] = FitFilter2Data(stim,response,OnlyThesePoints,flstr,regstr,regstr2);

		% make the prediction 
		fp = filter(K,1,stim-mean(stim)); % + mean(response);

		% apply the cutoff
		fp(fp<min_cutoff) = min_cutoff;

		% censor initial prediction 
		fp(1:filter_length+1) = NaN;

		% find the error--in r square
		err(i) = rsquare(fp(filter_length+2:end),response(filter_length+2:end));

		% find other metrics
		filter_height(i) = max(abs(K));
		filter_sum(i) = sum(abs(diff(K)));
		fall= fit(fp(filter_length+2:end),response(filter_length+2:end),'Poly1');

		mcond(i) = cond(C);
		slope(i) = fall.p1;
		

	end

	% save to diagnositics
	diagnostics.reg = reg;
	diagnostics.mcond = mcond;
	diagnostics.err = err;
	diagnostics.slope = slope;
	diagnostics.filter_sum = filter_sum;
	diagnostics.filter_height = filter_height;

	% picking the best filter
	if min(filter_sum) < filter_sum(1) && min(filter_sum) < filter_sum(end)
		% there is a minimum. just pick it
		[~,id]=min(filter_sum);
	else
		% pick a filter so that the error, the sum and slope are as close as possible to best values.
		err = abs(err - 1);
		filter_sum = (filter_sum - min(filter_sum))/min(filter_sum);
		slope = abs(1-slope);
		filter_sum = filter_sum/max(filter_sum);
		slope = slope/max(slope);
		[~,id]=min(max([filter_sum;err;slope]));
	end


	% and recalculate the filter
	flstr = strcat('filter_length=',mat2str(filter_length),';');
	regstr = strcat('reg=',mat2str(reg(id)),';');
	regstr2 = strcat('regtype=',mat2str(regtype),';');
	K = FitFilter2Data(stim,response,OnlyThesePoints,flstr,regstr,regstr2);
	diagnostics.bestfilter = id;

	% ensure unit gain
	K = K*diagnostics.slope(diagnostics.bestfilter);

else


	%% Damon's code
	nsteps = 20; 
	ss = (log(regmax)-log(regmin))/nsteps;
	reg = log(regmin):ss:log(regmax);
	reg = exp(reg);
	if length(reg) > nsteps
		reg(end) = [];
	end
	 

	err = NaN(1,nsteps);
	filter_height = NaN(1,nsteps);
	filter_sum = NaN(1,nsteps);
	slope = NaN(1,nsteps);

	for i = 1:nsteps
		% find the filter
		[Kdamon, ~, fp] = Back_out_1dfilter_new(stim,response,filter_length,reg(i));

		% correct the prediction 
		fp = fp+mean(response);

		% censor initial prediction 
		fp(1:filter_length+1) = NaN;

		% find the error
		err(i) = sqrt(sum((fp(filter_length+2:end) - response(filter_length+2:end)).^2 ));

		% find other metrics
		filter_height(i) = max(Kdamon);
		filter_sum(i) = sum(abs(Kdamon));
		[fall ~] = fit(fp(filter_length+2:end)',response(filter_length+2:end)','Poly1');
		slope(i) = fall.p1;

	end

	% log values to diagnostics
	diagnostics.reg = reg;
	diagnostics.err = err;
	diagnostics.slope = slope;
	diagnostics.filter_sum = filter_sum;
	diagnostics.filter_height = filter_height;

	% pick a filter so that the error, the sum and slope are as close as possible to best values.
	err = (err - min(err))/min(err);
	filter_sum = (filter_sum - min(filter_sum))/min(filter_sum);
	slope = abs(1-slope);
	[~,id]=min(max([filter_sum;err;slope]));

	% and recalculate the filter
	K = Back_out_1dfilter_new(stim,response,filter_length,reg(id));
	diagnostics.bestfilter = id;
end
