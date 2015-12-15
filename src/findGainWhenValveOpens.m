% findGainWhenValveOpens
% returns a gain vector evaluated at times when valve opens
% 
% created by Srinivas Gorur-Shandilya at 10:44 , 15 December 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function gain = findGainWhenValveOpens(valve,stim,resp)

assert(isvector(valve),'All arguments should be vectors')
assert(isvector(stim),'All arguments should be vectors')
assert(isvector(resp),'All arguments should be vectors')
assert(length(unique(valve))==2,'First argument should be a binary sequence')

% find valve openings
[ons,offs] = computeOnsOffs(valve);

% for each valve opening, find the peak
gain = NaN*stim;

for i = 1:length(ons)
	% find the maximum PID in the pulse
	temp = stim;
	temp(1:ons(i)-1) = 0;
	temp(offs(i):end) = 0;
	[max_pid,loc] = max(temp);
	gain(loc) = resp(loc)/max_pid;
end