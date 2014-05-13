% FitNonLinearity.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
% 
% fits a non-linear function, defined by a scaled hill function (with parameters A, K and n), to input xdata and output ydata
function [f A K n] = FitNonLinearity(xdata,ydata)
x0 = [10; 1; 2];
x = lsqcurvefit(@hill,x0,xdata,ydata);
A = x(1);
K = x(2);
n = x(3);
f = hill(x,xdata);

function [r] = hill(x,xdata)
	r = x(1)*(1./(1+(x(2)./xdata).^x(3)));