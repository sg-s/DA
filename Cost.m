% Cost.m
% computes a cost for two vectors
% created by Srinivas Gorur-Shandilya at 13:36 on 02-April-2014. Contact me at http://srinivas.gs/
function [c] = Cost(a,b)
a(1:2000) = [];
b(1:2000) = [];
a = a(:);
b = b(:);
c = sqrt(sum((a-b).^2));