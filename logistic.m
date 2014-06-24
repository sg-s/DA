
function [r] = logistic(x,xdata)
	A = x(1);
	b = x(2);
	k  =x(3);
	ee = b - k*xdata;
	r = A./(1 + exp(ee));
	