% created by Srinivas Gorur-Shandilya at 11:48 on 04-April-2014. Contact me at http://srinivas.gs/
% Dynamical Gain Model
% the Dynamical Gain model consists of two gamma filters, each characterised by theta and k, and a factor called p.f that specifies the mean baseline to be added to the output of the first gamma filter. 1 is added to the output of the second gamma filter, making this the "gain filter"
function [r] = dgm(x,p)
% unpack parameters
p1.theta = p.theta1;
p1.k = p.k1;
p2.theta = p.theta2;
p2.k = p.k2;
p3.theta = p.theta3;
p3.k = p.k3;

% process stimulus
x = x - mean(x);


x = x(:) - mean(x);

% make filters. all filters are a 1000 points long
t = 0:999;
K1 = GammaDist(t,p1) - GammaDist(t,p3);
K2 = GammaDist(t,p2);

% convolve stimulus with the gamma filters
y = filter(K1,1,x) + p.f;
g = filter(K1,1,x) + 1;

r = y.*g;