% FindBestFilter.m
% created by Srinivas Gorur-Shandilya at 17:31 , 15 January 2014. Contact me at http://srinivas.gs/contact/
% finds the best filter, searching through the space of the free parameter in my and damon's filter estimation functions
function [K Kdamon diagnostics] = FindBestFilter(stim,response,filter_length,CRange,DRange)

if nargin < 4
	regmax = 1e3;
	regmin = 1e-5;
 
else
	regmax = CRange(2);
	regmin = CRange(1);
end

%% Chichilnisky's method. 
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
	flstr = strcat('filter_length=',mat2str(filter_length),';');
	regstr = strcat('reg=',mat2str(reg(i)),';');
	K = FitFilter2Data(stim,response,flstr,regstr);

	% make the prediction 
	fp = filter(K,1,stim-mean(stim)) + mean(response);

	if size(fp,1) ~= 1
		fp = fp';
	end
	% censor initial prediction 
	fp(1:filter_length+1) = NaN;

	% find the error
	err(i) = sqrt(sum(((fp(filter_length+2:end)-response(filter_length+2:end)).^2)));

	% find other metrics
	filter_height(i) = max(K);
	filter_sum(i) = sum(abs(K));
	fall= fit(fp(filter_length+2:end)',response(filter_length+2:end)','Poly1');


	slope(i) = fall.p1;
	

end

% save to diagnositics
diagnostics.C.reg = reg;
diagnostics.C.err = err;
diagnostics.C.slope = slope;
diagnostics.C.filter_sum = filter_sum;
diagnostics.C.filter_height = filter_height;

% pick a filter so that the error, the sum and slope are as close as possible to best values.
err = (err - min(err))/min(err);
filter_sum = (filter_sum - min(filter_sum))/min(filter_sum);
slope = abs(1-slope);
[~,id]=min(max([filter_sum;err;slope]));

% and recalculate the filter
flstr = strcat('filter_length=',mat2str(filter_length),';');
regstr = strcat('reg=',mat2str(reg(id)),';');
K = FitFilter2Data(stim,response,flstr,regstr);
diagnostics.C.bestfilter = id;

if nargout == 1
	% no need to calculate Damon's filter
	return
end


if nargin < 5
	regmax = 1;
	regmin = 1e-7;
 
else
	regmax = CRange(2);
	regmin = CRange(1);
end


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
diagnostics.D.reg = reg;
diagnostics.D.err = err;
diagnostics.D.slope = slope;
diagnostics.D.filter_sum = filter_sum;
diagnostics.D.filter_height = filter_height;

% pick a filter so that the error, the sum and slope are as close as possible to best values.
err = (err - min(err))/min(err);
filter_sum = (filter_sum - min(filter_sum))/min(filter_sum);
slope = abs(1-slope);
[~,id]=min(max([filter_sum;err;slope]));

% and recalculate the filter
Kdamon = Back_out_1dfilter_new(stim,response,filter_length,reg(id));
diagnostics.D.bestfilter = id;
