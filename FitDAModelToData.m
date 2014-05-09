% FitDAModelToData.m
% fits the DA model to data
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
%
% data is a structure containing:
% .stimulus -- the stimulus
% .response -- the response of the ORN (or whatever sensor)
% .time -- a time vector
% 
% stimulus or response can be matrices, where each row represents a different experiment or a different trial. 
function [p Rguess] = FitDAModelToData(data)

x0 = [1200 16 0.3 0.66 3 6 2];
lb = [10 0 0 0 2 0 2];
ub = [5000 50 1 10 10 100 10];
psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',2000,'MaxFunEvals',10000);
x = patternsearch(@(x) DA_cost_function(x,data,@Cost2,'ga'),x0,[],[],[],[],lb,ub,psoptions);
p = ValidateDAParameters(x,'ga');


% debug
response = data.response;
stimulus = data.stimulus;
Rguess = 0*response;
for i = 1:size(stimulus,1)
	Rguess(i,:) = DA_integrate(stimulus(i,:),p);
end
clear i

Rguess = Rguess';
response = response';
Rguess =Rguess(:);
response = response(:);

time = data.time;

figure, multiplot([],response,Rguess)