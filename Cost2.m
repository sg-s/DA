% Cost.m
% computes a cost for two vectors
% created by Srinivas Gorur-Shandilya at 13:36 on 02-April-2014. Contact me at http://srinivas.gs/
function [c] = Cost2(a,b)

a = a(:);
b = b(:);


% penalise vectors that just aim for the mean or for very low values 

c = sqrt(sum((a-b).^2)); % distance to solution
