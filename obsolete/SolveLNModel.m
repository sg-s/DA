% SolveLNModel.m
% solves a Linear-Nonlinear model, i.e, a linear filter that operates on a static input non-linearity. 
% Usage:
% 1) generate output given a filter and a nonlinearity:
% [Rguess,K] = SolveLNModel(x,stimulus,R);
% where x is a 2-element vector specifying the shape of the hill input nonlinearity 
% the first element of x is the location of the inflection point
% and the second is the steepness
% stimulus is a long vector
% and R is a shorter vector that is the filter, in the same time units as stimulus. 
%
% if R is as long as stimulus, it is interpreted as the response, and a filter K is fit to the data. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [Rguess, K] = SolveLNModel(x,stimulus,R)	
if ~nargin
	help SolveNLModel
	return
end


if length(R) == length(stimulus)

	% R is the response; invert output function and get target output
	a = ihill2(x,R);

	% fit filter.
	K = FitFilter2Data(stimulus,a,[],'reg=0.1;','filter_length=201;');


else
	% use provided filter. 
	K = R;

end

a = filter(K,1,stimulus);

% throw out negative values
a(a<0)= 0;

% pass throught output non-linearity
Rguess = hill2(x,a);

