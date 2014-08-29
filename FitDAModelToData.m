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
% You can also specify which flavour of the DA model to use by specifying a global variable called DA_Model_Func which contains a function handle to the DA model you want to use. By default, DA_Model_Func is "@DA_integrate2"
%
% You also have to specify a corresponding Validate Param function handle in DA_Model_Validate_Param_Func
function [p, Rguess,x ] = FitDAModelToData(data,x0,lb,ub,IgnoreInitial)

global DA_Model_Func
global DA_Model_Validate_Param_Func	
global nsteps


switch nargin 
	case 0
		help FitDAModelToData
		return
	case 1
		error('Need more inputs. ')
	case 2
		lb = x0/2;
		ub = x0*10;
		IgnoreInitial = 300;
	case 3
		ub = x0*10;
		IgnoreInitial = 300;
	case 4
		IgnoreInitial = 300;
end


if isempty(nsteps)
	nsteps = 100;
end

psoptions = psoptimset('UseParallel',true, 'Vectorized', 'off','Cache','on','CompletePoll','on','Display','iter','MaxIter',nsteps,'MaxFunEvals',20000);

x = patternsearch(@(x) DA_cost_function(x,data,@Cost2,IgnoreInitial),x0,[],[],[],[],lb,ub,psoptions);


%psoptions = psoptimset('Display','iter','MaxIter',nsteps,'MaxFunEvals',20000);
%x = fminsearch(@(x) DA_cost_function(x,data,@Cost2,IgnoreInitial),x0,psoptions);


if isempty(DA_Model_Validate_Param_Func)
	DA_Model_Validate_Param_Func = @ValidateDAParameters2;
end

p = DA_Model_Validate_Param_Func(x);


stimulus = data.stimulus;
Rguess = 0*stimulus;

if isempty(DA_Model_Func)
	DA_Model_Func = @DA_integrate2;
end

if isvector(stimulus)
	Rguess = DA_Model_Func(stimulus,p);
else
	for i = 1:size(stimulus,2)
		Rguess(:,i) = DA_Model_Func(stimulus(:,i),p);
	end
	clear i
end

Rguess = abs(Rguess);