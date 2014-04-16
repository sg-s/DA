% Cost.m
% computes a cost for two vectors
% created by Srinivas Gorur-Shandilya at 13:36 on 02-April-2014. Contact me at http://srinivas.gs/
function [c] = Cost(a,b)

a = a(:);
b = b(:);

a = a(2000:end);
b = b(2000:end);

% penalise vectors that just aim for the mean or for very low values 

c1 = sqrt(sum((a-b).^2)); % distance to solution


m = mean(a); 
c2 = sqrt(sum((m-b).^2));  % distance to mean

c = c1/c2;