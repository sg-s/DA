% findValidHistoryLengths.m
% finds valid history lengths based on stimulus, prediction and response so that green and red dots have a good degree of overlap
% 
% created by Srinivas Gorur-Shandilya at 9:59 , 16 July 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [history_lengths] = findValidHistoryLengths(dt,stim,pred,resp,n_points,frac);

overlap_fraction = .1;

% prepare the data
x.stimulus = stim;
x.frac = frac;
x.response = resp;
x.prediction = pred;
x.time = dt*(1:length(stim));

% look only in these history lengths
history_lengths = logspace(-1,0,30);
example_history_length = history_lengths(10);

p = NaN(2,length(history_lengths));

[p,low_slopes,high_slopes,low_gof,high_gof,example_plot,extra_variables] = GainAnalysis5(x,history_lengths,example_history_length,[],p);

% find the overlap
ok1 = (extra_variables.low_max - extra_variables.high_min)./(extra_variables.low_max-extra_variables.low_min)>overlap_fraction;
 ok2 = (extra_variables.low_max - extra_variables.high_min)./(extra_variables.high_max-extra_variables.high_min)>overlap_fraction;
 ok = ok1 & ok2;

 min_h = history_lengths(find(ok,1,'first'));

 max_h = dt*(length(stim)/5);

 history_lengths = logspace(log10(min_h),log10(max_h),n_points);

