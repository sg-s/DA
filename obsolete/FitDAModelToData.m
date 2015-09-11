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
function [p, Rguess,x ] = FitDAModelToData(data,x0,lb,ub,min_r2)

default_x0 = [max(data.response) std(data.response) .3 10 2 10 2 -.1*(max(data.stimulus))];
scale = 4;

switch nargin 
	case 0
		help FitDAModelToData
		return
	case 1
		x0 = default_x0;
		lb = x0/scale;
		ub = x0*scale;
		min_r2 = .95;
	case 2
		if isempty(x0)
			x0 = default_x0;
		end
		lb = x0/scale;
		ub = x0*scale;
		min_r2 = .95;
	case 3
		ub = x0*scale;
	case 4
	case 5
		if isempty(x0)
			x0 = default_x0;
		end
		if isempty(lb)
			lb = x0/scale;

		end
		if isempty(ub)
			ub = x0*scale;
		end

end

% special bounds
temp=ub(x0<0);
ub(x0<0) = lb(x0<0);
lb(x0<0) = temp; clear temp;
lb(3) = 0; ub(3) = 1;

IgnoreInitial = 300;
nrep = 20;


nsteps = 70;


psoptions = psoptimset('UseParallel',true, 'Vectorized', 'off','Cache','on','CompletePoll','on','Display','iter','MaxIter',nsteps,'MaxFunEvals',20000);


DA_Model_Validate_Param_Func = @ValidateDAParameters2;
DA_Model_Func = @DA_integrate2;

stimulus = data.stimulus;
Rguess = 0*stimulus;


if min_r2
	psoptions = psoptimset('UseParallel',true, 'Vectorized', 'off','Cache','on','CompletePoll','on','Display','iter','MaxIter',nsteps,'MaxFunEvals',20000);
	x = patternsearch(@(x) DA_cost_function(x,data,IgnoreInitial),x0,[],[],[],[],lb,ub,psoptions);

	% keep crunching till we can fit the damn thing
	x = x0;
	
	for i = 1:nrep

		% specify new bounds
		lb = x/scale;
		ub = x*scale;
		% special bounds
		temp=ub(x<0);
		ub(x<0) = lb(x<0);
		lb(x<0) = temp; clear temp;
		if lb(3) < 0
			lb(3) = 0; 
		end
		if ub(3) > 1
			ub(3) = 1;
		end
		
		% make sure that p.s0 does not exceed the minimum value of the stimulus
		if ub(end) > min(stimulus)
			ub(end) = min(stimulus)-1e-5;
		end
		if lb(end) > ub(end)
			lb(end) = 0;
		end
		if x(end) > ub(end)
			x(end) = ub(end)/2;
		end

		% fit
		x = patternsearch(@(x) DA_cost_function(x,data,IgnoreInitial),x,[],[],[],[],lb,ub,psoptions);
		p = DA_Model_Validate_Param_Func(x);
		
		if isvector(stimulus)
			Rguess = DA_Model_Func(stimulus,p);
		else
			for i = 1:size(stimulus,2)
				Rguess(:,i) = DA_Model_Func(stimulus(:,i),p);
			end
			clear i
		end

		Rguess = abs(Rguess);

		if rsquare(Rguess(IgnoreInitial:end),data.response(IgnoreInitial:end)) > min_r2
			return
		else
			disp(oval(rsquare(Rguess(IgnoreInitial:end),data.response(IgnoreInitial:end)),4))
			fprintf('\n')
		end	

	end
else
	x = patternsearch(@(x) DA_cost_function(x,data,IgnoreInitial),x0,[],[],[],[],lb,ub,psoptions);

end


p = DA_Model_Validate_Param_Func(x);


if isvector(stimulus)
	Rguess = DA_Model_Func(stimulus,p);
else
	for i = 1:size(stimulus,2)
		Rguess(:,i) = DA_Model_Func(stimulus(:,i),p);
	end
	clear i
end

Rguess = abs(Rguess);

% show the result
figure, hold on
plot(data.response,'k')
plot(Rguess,'r')
title(oval(rsquare(Rguess(IgnoreInitial:end),data.response(IgnoreInitial:end)),4))