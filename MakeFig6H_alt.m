% MakeFig6H_alt.m
% an alternative version of Make Fig 6H, where we rely on peaks in the signal and response
% to measure temporal changes in dynamics. 
% 
% created by Srinivas Gorur-Shandilya at 1:21 , 05 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function c = MakeFig6H_alt(stimulus,response,history_length)
if ~nargin
	help MakeFig6H
	return
end
if ~isvector(stimulus) || ~isvector(response) 
	error('First two arguments should be vectors')
else
	stimulus = stimulus(:);
	response =  response(:);
	if length(stimulus) ~= length(response)
		error('Stimulus and response have to have the same length')
	end
end
if length(history_length) > 1
	error('History Length has to be a scalar')
end

c = NaN*stimulus;

stimulus = stimulus - min(stimulus);
stimulus = stimulus/max(stimulus);
response = response - min(response);
response = response/max(response);

[~,loc] = findpeaks(stimulus,'MinPeakProminence',.002);

for i = 1:length(loc)
	snippet = response(loc(i):loc(i)+300);
	[~,c(loc(i))]=max(snippet);
end

c(c>300) = NaN;