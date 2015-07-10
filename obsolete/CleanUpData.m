% CleanUpData.m
% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
% cleans up data (as specified in FitDAModelToData.m)
% so that:
% stimulus is positive, and is zero when it should be
function [data] = CleanUpData(data)

response = data.response;
stimulus = data.stimulus;
time = data.time;

% signal starts here
stimulusON = 1;

sz = size(stimulus);
for i = 1:sz(1)
	a = find(time<stimulusON,1,'last');
	offset = mean(stimulus(i,1:a));
	stimulus(i,:) = stimulus(i,:) - offset;
end
clear i

stimulus = abs(stimulus);

% return
data.stimulus = stimulus;
