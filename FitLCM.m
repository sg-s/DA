% FitLCM.m
% Fits a Linearised Cascade Model to data
%
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [] = FitLCM(data,model,max_order,x0)
switch nargin
case 0
	help FitLCM
	return
case 1
	model = @gamma_function;
	max_order = 2;
	x0 = [1 1 1];
case 2
	% assume that the model requires N parameters, and one time input
	n = nargin(model);
	x0 = ones(1,n);
case 3
end

% parameters
nsteps = 500;
IgnoreInitial = 201;

% build a simple LN model 
td=1;
[K,~,filtertime] = FindBestFilter(data(td).PID(500:end),data(td).ORN(500:end),[],'filter_length=201;');
data(td).K = K;
data(td).filtertime = filtertime*mean(diff(data(td).time));
data(td).LinearFit = convolve(data(td).time,data(td).PID,data(td).K,data(td).filtertime);
data(td).LinearFit = data(td).LinearFit + mean(data(td).ORN);

xdata = data(td).LinearFit;
ydata = data(td).ORN;

% crop it to lose NaNs
ydata(isnan(xdata)) = [];
xdata(isnan(xdata)) = [];

xdata = xdata(:);
ydata = ydata(:);

fo=optimset('MaxFunEvals',1000,'Display','none');
x = lsqcurvefit(@hill,[max(ydata) 2 2],xdata,ydata,[max(ydata)/2 2 1],[2*max(ydata) max(ydata) 10],fo);


% save this for later
LNFit = hill(x,data(td).LinearFit);



psoptions = psoptimset('UseParallel',true, 'Vectorized', 'off','Cache','on','CompletePoll','on','Display','final','MaxIter',nsteps,'MaxFunEvals',20000);


% debug
colours = jet(max_order+1);
figure, plot(data.ORN(IgnoreInitial:end),'k')
hold on
plot(LNFit(IgnoreInitial:end),'Color',colours(1,:))


% now fit the cascade model
rcap = NaN(max_order+1,length(data.ORN(IgnoreInitial:end)));
rcap(1,:) = LNFit(IgnoreInitial:end); % stores the nth order predictions
fit_cost = NaN(1,length(max_order)+1);
fit_param = NaN(length(x0),length(max_order));


numerator = LNFit(IgnoreInitial:end);
x0 = [10 .2 10];
d.stimulus = data.PID-mean(data.PID);
d.stimulus = d.stimulus(IgnoreInitial:end);
d.response = data.ORN(IgnoreInitial:end);


for i = 1:max_order
	% assemble the data
	d.numerator = numerator;
	ub = x0*100; lb = x0/100;

	% fit
	x = patternsearch(@(x) LCM_Cost_Function(x,d,model),x0,[],[],[],[],lb,ub,psoptions);
	% save for later
	disp(x)
	fit_param(:,i) = x;

	% debug
	K = model(x(1),x(2),x(3));
	Rguess = filter(K,1,d.stimulus);
	Rguess = d.numerator./(1+Rguess);
	rcap(1+i,:) = Rguess;
	plot(Rguess,'Color',colours(1+i,:))

	numerator = Rguess;
	

end

% calculate the cost at each iteration
 a = data.ORN(IgnoreInitial:end);
for i = 1:width(rcap)
	fit_cost(i) = Cost2(a(~isnan(rcap(i,:))),rcap(i,(~isnan(rcap(i,:)))));
end


% make a figure showing all the gain filters calcualted
figure, hold on
for i = 1:max_order
	K = model(fit_param(:,i));
	plot(K,'Color',colours(i+1,:))
end

keyboard







