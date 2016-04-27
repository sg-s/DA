% computeSmoothedStimulus
% part of the DA project
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [shat] = computeSmoothedStimulus(stimulus,hl)
switch nargin 
	case 0
		help computeSmoothedStimulus
		return
	case 1
		error('Not enough input arguments')
	case 2
end

assert(~any(isnan(hl)), 'computeSmoothedStimulus: History lengths are NaN.')
assert(~any(isnan(stimulus)), 'computeSmoothedStimulus: Stimulus has NaNs.')


shat = NaN(length(hl),length(stimulus));
for i = 1:length(hl)
	if hl(i) == 0
		shat(i,:) = stimulus;
	else
		shat(i,:) = filter(ones(1,hl(i))/hl(i),1,stimulus);
		shat(i,1:hl(i)) = NaN;
	end
 	
end
