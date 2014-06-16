


filename = '/data/orn/random-stim/final_2011_06_03_ab3A_1but5X-4_20ml_30sec_100ms_rand.mat';
[PID, time, f,Valve,uncropped] = PrepData3(filename);
PID = PID(:);
time = time(:);
f = f(:);

% detrend PID
ptrend = fit(time,PID,'Poly1'); 
PID = PID - (ptrend(time) - mean(ptrend(time)));


% assemble into a data structure
data(1).stimulus = Valve;
data(1).response = f;
data(1).Valve = Valve;
data(1).time = time;


% build a simple linear model
[K,~,filtertime] = FindBestFilter(PID(500:end),f(500:end),[],'filter_length=201;');
data(1).K = K;
LinearFit = filter(K,1,PID-mean(PID))+ mean(f(500:end));
LinearFit(LinearFit<min(f)) = min(f);
% correct for offset
offset = abs(min(filtertime));
LinearFit(1:offset) = [];
LinearFit = [LinearFit; NaN(offset,1)];

% add it to the data
data(1).LinearFit = LinearFit;

figure('outerposition',[0 0 1500 1000],'PaperUnits','points','PaperSize',[1000 500]); hold on
fig(1).a(1) = subplot(3,6,2:6);
fig(1).a(2) = subplot(3,6,8:12);
fig(1).a(3) = subplot(3,6,13);
fig(1).a(4) = subplot(3,6,14:18); hold on

filtertime = filtertime*mean(diff(time));
offset = offset*mean(diff(time));

plot(fig(1).a(1),time,Valve);
plot(fig(1).a(2),time,PID);
plot(fig(1).a(3),filtertime,K);
plot(fig(1).a(4),time,f,'k'); 
plot(fig(1).a(4),time,LinearFit,'r')
set(fig(1).a(4),'XLim',[20 25])
set(fig(1).a(1),'XLim',[20 25],'YLim',[-0.1 1.1])
set(fig(1).a(2),'XLim',[20 25])
set(fig(1).a(3),'XLim',[min(filtertime-offset) max(filtertime)+offset])
PrettyFig;


return


	x0 = [27596       572  0.1  5   2   15  27  0];
	lb = [100         0    0     0   2   0   2  0];
	ub = [89200       1900  1    10  2  30   2  10];
	psoptions = psoptimset('UseParallel',true,'CompletePoll', 'on', 'Vectorized', 'off','Display','iter','MaxIter',2000,'MaxFunEvals',50000);
	x = patternsearch(@(x) DA_cost_function(x,data(1),@Cost,500,0),x0,[],[],[],[],lb,ub,psoptions);
	data(1).DAFit_p = ValidateDAParameters2(x);
	pred = DA_integrate2(data(1).stimulus,data(1).DAFit_p);
	multiplot([],data(1).response,pred)
	clear x
	data(1).DAFit_p
	% save('DA_Paper_data.mat','data','-append')

