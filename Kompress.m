% created by Srinivas Gorur-Shandilya at 13:54 on 02-April-2014. Contact me at http://srinivas.gs/
% Kompress accepts any number, and passes it through a function whose range is [0,1]
% Kompress can also invert the function, if called with 'invert'
function [k] = Kompress(a,action)
if nargin == 1
	action = 'Kompress';
end
switch action 
case 'Kompress'
	k = a/(sqrt(1+a^2));
	k = (1+k)/2;
case 'invert'
	error
end