% created by Srinivas Gorur-Shandilya at 11:48 on 04-April-2014. Contact me at http://srinivas.gs/
% Dynamical Gain Model
function [r] = dgm(time,x,p)

x = x(:) - mean(x);
dt = mean(diff(time));

% make filters
K1 = 

% convolve stimulus with first gamma filter
y = filter();