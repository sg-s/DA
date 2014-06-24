% FitNonLinearity.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
% 
% fits a non-linear function, defined by a scaled logistic function (with parameters A, K and n), to input xdata and output ydata
function [f A K n] = FitNonLinearity(xdata,ydata,model)

yl = length(ydata);

if nargin < 3
	model = 'logistic';
end

xdata = xdata(:);
ydata = ydata(:);

% weed out NaNs
rm1 = isnan(xdata);
rm2 = isnan(ydata);
rm_this = logical(rm1+rm2);


xdata(rm_this) = [];
ydata(rm_this) = [];

x0 = [10; 1; 2];
options = optimoptions('lsqcurvefit','MaxFunEvals',10000,'MaxIter',1000);
switch model
	case 'hill'
		x = lsqcurvefit(@hill,x0,xdata,ydata,[0 0 1],[],options);
		f = hill(x,xdata);
	case 'logistic'
		x = lsqcurvefit(@logistic,x0,xdata,ydata,[],[],options);
		f = logistic(x,xdata);
end

A = x(1);
K = x(2);
n = x(3);

% pad back with NaNs
if length(f) < yl
	pl = yl - length(f);
	f = [f(:); NaN(pl,1)];
end

