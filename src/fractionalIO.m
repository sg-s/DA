% fractionalIO.m
% makes a fractional I/O curve, to estimate gain in a non-dimensional way
% 
% created by Srinivas Gorur-Shandilya at 5:40 , 03 October 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [x,y] = fractionalIO(time,stim,resp,filtertime,K)

% defensive programming 
assert(isvector(stim) & isvector(resp), 'First two arguments should be vectors')
assert(length(stim) == length(resp), 'First two argument should have equal length')
assert(isvector(K),'Third argument should be a filter, a vector')

stim = stim(:);
resp = resp(:);
K = K(:);

x = convolve(time,stim,K,filtertime) + nanmean(stim);

x = (x - nanmean(x))./nanmean(x);
y = (resp - nanmean(resp))./nanmean(resp);
