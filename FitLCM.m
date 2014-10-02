% FitLCM.m
% Fits a Linearised Cascade Model to data
%
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [data] = FitLCM(data,model,max_order,x0)

param_max_order = 4;	
switch nargin
case 0
	help FitLCM
	return
case 1
	model = @alpha_filter;
	max_order = param_max_order; % max order of gain correction
	x0 = [10 1];
case 2
	
	max_order = param_max_order;
case 3
	n = nargin(model)-1;
	x0 = rand(1,n);
end

if nargin < 4
	switch char(model)
		case 'filter_gamma'
			x0 = [10 2 1];
			ub = [100 4 1e5];
			lb = [1e-3 0 0];
		case 'filter_alpha'
			x0 = [10 1];
			ub = [1e2 1e4];
			lb = [1 1e-4];
		case 'filter_alpha2'
			x0 = [30 60 1 .1];
			ub = [300 300 1e4 2  ];
			lb = [1/10 1/10 0 ];
		case 'filter_gamma2'
			x0 = [10 2 200 1     ];
			ub = [1e3 5 1e3 1e4  ];
			lb = [1 1 2 1        ];
	end
end

% parameters
nsteps = 2000;
IgnoreInitial = 201;
IgnoreTerminal = 102;

% build a simple LN model if not provided in data
if isfield(data,'f0')
else
	[K,~,filtertime] = FindBestFilter(data.PID(500:end),data.ORN(500:end),[],'filter_length=201;');
	data.K = K;
	data.filtertime = filtertime*mean(diff(data.time));
	data.LinearFit = convolve(data.time,data.PID,data.K,data.filtertime);
	data.LinearFit = data.LinearFit + mean(data.ORN);

	xdata = data.LinearFit;
	ydata = data.ORN;

	% crop it to lose NaNs
	ydata(isnan(xdata)) = [];
	xdata(isnan(xdata)) = [];

	xdata = xdata(:);
	ydata = ydata(:);

	fo=optimset('MaxFunEvals',1000,'Display','none');
	x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);


	% save this for later
	data.f0 = hill(x,data.LinearFit);
end


psoptions = psoptimset('UseParallel',true, 'Vectorized', 'off','Cache','on','CompletePoll','on','Display','final','MaxIter',nsteps,'MaxFunEvals',20000);


% debug
if ~nargout
	colours = jet(max_order+1);
	figure, plot(data.ORN(IgnoreInitial:end-IgnoreTerminal),'k')
	hold on
	plot(data.f0(IgnoreInitial:end-IgnoreTerminal),'Color',colours(1,:))
end


% now fit the cascade model
rcap = NaN(max_order+1,length(data.ORN(IgnoreInitial:end-IgnoreTerminal)));
rcap(1,:) = data.f0(IgnoreInitial:end-IgnoreTerminal); % stores the 0th order predictions
fit_cost = NaN(1,length(max_order)+1);
fit_param = NaN(length(x0),length(max_order));


numerator = data.f0(IgnoreInitial:end-IgnoreTerminal);
d.stimulus = data.PID;
d.stimulus = d.stimulus(IgnoreInitial:end-IgnoreTerminal);
d.stimulus = d.stimulus - mean(d.stimulus);
d.response = data.ORN(IgnoreInitial:end-IgnoreTerminal);


for i = 1:max_order
	% assemble the data
	d.numerator = numerator;

	% fit
	x = patternsearch(@(x) LCM_Cost_Function(x,d,model),x0,[],[],[],[],lb,ub,psoptions);
	% save for later
	disp(x)
	fit_param(:,i) = x;

	% debug

	n = nargin(model)-1;
	K = model(x(1:n));
	Rguess = filter(K,1,d.stimulus);
	Rguess = d.numerator./(1+Rguess);
	rcap(1+i,:) = Rguess;
	if ~nargout
		plot(Rguess,'Color',colours(1+i,:))
	end

	numerator = Rguess;


end

keyboard
% calculate the cost at each iteration
a = data.ORN(IgnoreInitial:end-IgnoreTerminal);
for i = 1:width(rcap)
	fit_cost(i) = Cost2(a(~isnan(rcap(i,:))),rcap(i,(~isnan(rcap(i,:)))));
end


% make a figure showing all the gain filters calcualted
if ~nargout
	figure, hold on
	for i = 1:max_order
		n = nargin(model)-1;
		K = model(fit_param(1:n,i));
		plot(K,'Color',colours(i+1,:))
	end


	% make a figure showing the improvement in fit with order
	figure, hold on
	plot(fit_cost,'k')


	tilefigs;
end





