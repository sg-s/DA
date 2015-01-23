% SubtractScale.m
% this subtracts a number and scales a vector by another to fit some target
% a very simple model indeed
%
% created by Srinivas Gorur-Shandilya at 3:31 , 23 January 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function f = SubtractScale(s,p)
s = s + p.A;
f = s*p.B;