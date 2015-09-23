% FitLinearModel2Data.m
% accepts a filter, and a dataset, and finds the best scaling of the filter to match the data. 
% data is a structure with fields called stimulus and response and time. 
% make sure the filter has the same time units as the data. 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [LinearFit Kscaled] = FitLinearModel2Data(data,K)




x0 = [1/100];
lb = [1e-10];
ub = [10];

psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',300,'MaxFunEvals',1000);
x = patternsearch(@(x) LinearModelCostFunction(x,data,K),x0,[],[],[],[],lb,ub,psoptions);

% unpack data
stimulus = data.stimulus(:);
response = data.response(:);

% offset
offset = floor(length(K)/10);
stimulus = stimulus(offset:end);
response = response(1:end-offset+1);


%scale filter
K = K*x;

% filter
Rguess = (filter(K,1,stimulus-mean(stimulus))+mean(response));
Rguess(Rguess<0)=0;
% disp(x)
% figure, hold on
% plot(response)
% plot(Rguess,'r')

LinearFit = [Rguess; zeros(offset-1,1)];
Kscaled = K;


function [c] = LinearModelCostFunction(x,data,K);

% unpack data
stimulus = data.stimulus(:);
response = data.response(:);

% offset
offset = floor(length(K)/10);
stimulus = stimulus(offset:end);
response = response(1:end-offset+1);

%scale filter
K = K*x;

% filter
Rguess = (filter(K,1,stimulus-mean(stimulus))+mean(response));

Rguess(Rguess<0)=0;


% compute cost
c = Cost(Rguess,response);

