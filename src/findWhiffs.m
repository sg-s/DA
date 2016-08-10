% findWhiffs.m
% find whiffs in stimulus as defined by regions when the stimulus sharply increases
% 
% created by Srinivas Gorur-Shandilya at 1:30 , 15 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [whiff_ons,whiff_offs] = findWhiffs(s)


min_whiff_duration = 10; % time units, probably ms

s = filter(ones(10,1),1,s);
ds = diff(s); ds(ds<0) = 0; ds=ds/std(ds); ds(ds>0) = 1;
[whiff_ons,whiff_offs] = computeOnsOffs(ds);

rm_this = (whiff_offs-whiff_ons) < min_whiff_duration;
whiff_ons(rm_this) = [];
whiff_offs(rm_this) = [];