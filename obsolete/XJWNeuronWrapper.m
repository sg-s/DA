% implements the Liu-Wang model of spike-based adaptation
% 
% created by Srinivas Gorur-Shandilya at 9:24 , 23 April 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [V,f,Ca] = XJWNeuronWrapper(time,stimulus,p)

stimulus = stimulus(:);
time = time(:);

nrep = 10;
noise = 0;


st = zeros(length(time),nrep);


for i = 1:nrep
	[V,Ca,st(:,i)] = XJWNeuronEuler(time,stimulus+noise*randn(length(time),1),p);
end


f=spiketimes2f(st,time,1e-3,0.03);
t = min(time):1e-3:max(time);


if ~nargout
	figure, hold on
	plot(f)
else
	% interpolate to get back to the original length
	f = mean(f,2);
	f = interp1(t,f,time);
end 