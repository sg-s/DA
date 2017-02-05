%%
%  In this document, I fit a single adapting NLN model (aNLN4) to all the data.

% get the filter from the Gaussian stimuli 
clear MSGdata
MSGdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
MSGdata = cleanMSGdata(MSGdata);

% make sure stimulus is always positive
for i = 1:size(MSGdata.PID,2)
	MSGdata.PID(:,i) = MSGdata.PID(:,i) - min(MSGdata.PID(:,i));
end

a = 25e3; z = 55e3;
clear data
for i = 1:max(MSGdata.paradigm)
	S = MSGdata.PID(:,MSGdata.paradigm==i);
	R = MSGdata.fA(:,MSGdata.paradigm==i);
	rm_this = sum(R) == 0;
	S(:,rm_this) = []; R(:,rm_this) = [];
	data(i).stimulus = mean(S(a:z,:),2);
	data(i).response = mean(R(a:z,:),2);
	data(i).response(1:10e3) = NaN;
end


clear p
p.    k0 = 0.0968;
p.     B = 0.8289;
p.tau_y1 = 27.1055;
p.tau_y2 = 96.9453;
p.     A = 0.7000;
p.     C = 178.9440;
p.     n = 2;


figure('outerposition',[0 0 1101 901],'PaperUnits','points','PaperSize',[1101 901]); hold on
clear ax
ax(1) = subplot(4,3,1); hold on
ax(2) = subplot(4,3,2); hold on
ax(3) = subplot(4,3,3); hold on


% first, show 3 nonlinearities from the Gaussians with the right colours 
c = parula(11);
for i = [2 6 10]
	this_paradigm = find(MSGdata.paradigm == i,1,'first');
	S = MSGdata.PID(:,this_paradigm);
	[R, a, b, kD, K] = aNLN4(S,p);
	% plot actual curve with mean kD
	x = logspace(-2,2,100);
	A = 1./(1 + mean(kD(35e3:55e3))./x);
	plot(ax(1),x,A,'Color',c(i,:))
	% now plot the points on it that are actually achieved 
	plot(ax(1),S(35e3:55e3),a(35e3:55e3),'.','Color',c(i,:))
end
set(ax(1),'XScale','log','XLim',[1e-2 10])
ax(1).YAxisLocation = 'right';
xlabel(ax(1),'Stimulus (V)')
ylabel(ax(1),'a')
axes(ax(1))
axis square

% show the filter 
plot(ax(2),K/norm(K),'k')
set(ax(2),'XLim',[0 600])
xlabel(ax(2),'Lag (ms)')
ylabel(ax(2),'Filter (norm)')
axes(ax(2))
axis square

% show the naturalistic stimulus fits

cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);
time = 1e-3*(1:length(data(1).S));

% show the time series of the natualistic stimulus
ax(4) = subplot(4,3,4:5); hold on
clear l
l(1) = plot(ax(4),time,data(2).R(:,3),'k');
% generate responses and plot them too

for j = 1:size(data(2).X,2)
	data(2).P(:,j) = aNLN4(data(2).S(:,j) - min(data(2).S(:,j)) ,p);
end
% fix some trivial scaling
ff = fit(data(2).P(:),data(2).R(:),'poly1');
for j = 1:size(data(2).X,2)
	data(2).P(:,j) = ff(data(2).P(:,j));
end

l(2) = plot(ax(4),time,data(2).P(:,3),'r');
legend(l,{'ab2A',['model r^2 = ' oval(rsquare(data(2).P(:,3),data(2).R(:,3)))]})

xlabel(ax(4),'Time (s)')
ylabel(ax(4),'Firing rate (Hz)')



