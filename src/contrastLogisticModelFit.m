% contrastLogisticModelFit
% this model is a wrapper around contrastLogisticModel
% used to fit the model
% this fits directly to the input output curves of a contrast-switch experiment
% so the stimulus construction is very weird and convoluted.
% S is a matrix such that:
% the 2nd-nth rows contain many trials of the stimulus,
% where
% each trial is evenly split between first a hi variance trial and then a lo variance trial
% and 
% the first row has the dose-response info as follows:
% [hi_x, lo_x] and is padded with NaNs. 
% 
% created by Srinivas Gorur-Shandilya at 1:49 , 16 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [R] = contrastLogisticModelFit(S,p)

% for logistic function
p.x0;
p.A;

% for changing k:
p.k0;
p.n;
p.tau;
p.B;

% bounds
lb.tau = 1;
lb.n = 1;
lb.B = eps;
lb.A = 57;
lb.k0 = 0;
lb.x0 = -1.325;
ub.x0 = -1.325;

ub.tau = 200;
ub.n = 4;
ub.B = 1000;

% strip out the dose-response info
temp = S(:,1);
S = S(:,2:end);

temp = nonnans(temp);
w = round(length(temp)/2);
hi_x = temp(1:w);
lo_x = temp(w+1:end);

% for each trial, compute the steepness of the contrast logisitic model
k = NaN*S;

for i = 1:width(S)
	[~,~,k(:,i)] = contrastLogisticModel(S(:,i),p);
end

% assumes 10-second long trials
k_hi = nanmean(nanmean(k(1e3:5e3,:)));
k_lo = nanmean(nanmean(k(6e3:end,:)));


y_hi = logistic(hi_x,p.A,k_hi,p.x0);
y_lo = logistic(lo_x,p.A,k_lo,p.x0);

R = [y_hi; y_lo];