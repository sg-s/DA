% function [rt] = FindResponseTimes(x)
% finds the response times (in matrix indices of a given matrix x, after specifing when the stimulus nominally starts with t, also a matrix index)
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
% 
function [rt] = FindPeakResponseTimes(x,t)
% x should be time by ntrials
rt = NaN(size(x,2),1);
m = max(x);
for i = 1:size(x,2)
	rt(i)= find(x(:,i)==m(i),1,'first') - t;
end
clear i

% rt can't be negative
rt(rt<0)=NaN;