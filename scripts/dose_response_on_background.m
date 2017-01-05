
%% Dose Responses on backgrounds
% In this document, I plot responses of the ab3 sensillum to ethyl acetate pulses, either on no background, or on top of backgrounds of ethyl acetate. The point is to directly observe the change in the dose-response curve of the neuron with adaptation. 

pHeader;

% gather data

root = '/data/DA-paper/data-for-paper/dose-response-on-bg/'; 
allfiles = dir([root '*.kontroller']);
allfiles(cellfun(@(x) strcmp(x(1),'.'), {allfiles.name})) = [];

S = [];
R = [];
fly = [];
sensillum = [];
back_level = [];

for i = 1:length(allfiles)
	load([root allfiles(i).name],'-mat')

	this_fly = str2double(allfiles(i).name(strfind(allfiles(i).name,'_F')+2:strfind(allfiles(i).name,'_F')+2));
	this_back_level = str2double(allfiles(i).name(strfind(allfiles(i).name,'_b')+2:strfind(allfiles(i).name,'_b')+2));
	this_sensillum = str2double(allfiles(i).name(strfind(allfiles(i).name,'_S')+2:strfind(allfiles(i).name,'_S')+2));
	this_S = vertcat(data.PID)'; this_S = this_S(1:10:end,:);
	this_R = vertcat(data.voltage)'; this_R = this_R(1:10:end,:);

	this_fly = this_fly*ones(size(this_S,2),1);
	this_back_level = this_back_level*ones(size(this_S,2),1);
	this_sensillum = this_sensillum*ones(size(this_S,2),1);

	% consolidate
	fly = [fly; this_fly];
	sensillum = [sensillum; this_sensillum];
	back_level = [back_level; this_back_level];
	S = [S this_S];
	R = [R this_R];
end

% for each trial, compute metrics of interest 
dR = NaN*this_back_level;
dS = NaN*this_back_level;
R_baseline = NaN*this_back_level;
S_baseline = NaN*this_back_level;

% correct stimulus to remove the baseline
all_orns = unique(sensillum);
for i = 1:length(all_orns)
	min_S = min(min(S(1:1e3,sensillum == all_orns(i) & back_level == 0)));
	S(:,sensillum == all_orns(i)) = S(:,sensillum == all_orns(i)) - min_S;
end

for i = 1:size(R,2)
	R_baseline(i) = mean(R(1:1e3,i));
	dR_plus = max(R(1e3:1.5e3,i)) - R_baseline(i);
	dR_minus = min(R(1e3:1.5e3,i)) - R_baseline(i);
	if abs(dR_plus) > abs(dR_minus)
		dR(i) = dR_plus;
	else
		dR(i) = dR_minus;
	end
	S_baseline(i) = mean(S(1:1e3,i));
	dS(i) = mean(S(1.1e3:1.5e3,i)) - S_baseline(i);
end

% ft = fittype('-(A.*((x).^n)./((x).^n + K.^n) + B)');
% ft2 = fittype('-(A.*(exp(x).^n)./(exp(x).^n + exp(K).^n) + B)');

%%
% In the following two figure, I plot (a) the statistics of the background stimulus (before foreground pulse presentation).  Note that there is some variability in the background level due to errors in odor delivery. (b) Change in LFP response (peak - baseline) vs. change in stimulus (peak - baseline). (c) Absolute LFP response vs. stimulus. (d) Change in response vs. peak stimulus. In all plots,points are coloured and grouped by target background level. Each figure comes from a single sensillum. 


for i = 1:length(all_orns)
	figure('outerposition',[0 0 900 905],'PaperUnits','points','PaperSize',[900 905]); hold on
	subplot(2,2,1); hold on
	all_back_levels = unique(back_level(sensillum == all_orns(i)));
	c = parula(length(all_back_levels)+1);
	for j = 1:length(all_back_levels)
		plot_these = sensillum == all_orns(i) & back_level == all_back_levels(j);
		plot(mean(S(1:1e3,plot_these)),std(S(1:1e3,plot_these)),'.','Color',c(j,:),'MarkerSize',20)
	end
	xlabel('\mu_{stimulus} (V)')
	ylabel('\sigma_{stimulus} (V)')
	title('Background stimulus statistics')

	% plot change in response vs. change in stimulus
	subplot(2,2,2); hold on
	for j = 1:length(all_back_levels)
		plot_these = sensillum == all_orns(i) & back_level == all_back_levels(j);
		x = dS(plot_these);
		y = dR(plot_these);
		[x,idx] = sort(x);
		plot(x,y(idx),'.-','Color',c(j,:),'MarkerSize',20)
		y(x<0) = []; x(x<0) = [];
		% if length(x) > 4
		% 	%ff = fit(log(x),y,ft2,'StartPoint',[1 0 mean((x)) 1],'Lower',[0 -10 -Inf 0],'Upper',[5 5 2*max(x) 2]);
		% 	ff = fit(x,y,ft,'StartPoint',[-1 0],'Upper',[ 0 0]);
		% 	plot(logspace(-4,2,100),ff(logspace(-4,1,100)),'Color',c(j,:))
		% end
	end
	xlabel('\DeltaS')
	ylabel('\DeltaR')
	set(gca,'XScale','log','XLim',[1e-4 10],'XTick',[1e-4 1e-3 1e-2 1e-1 1 10])

	% plot the response vs. stimulus
	subplot(2,2,3); hold on
	for j = 1:length(all_back_levels)
		plot_these = sensillum == all_orns(i) & back_level == all_back_levels(j);
		x = dS(plot_these)+S_baseline(plot_these);
		y = dR(plot_these)+R_baseline(plot_these);
		[x,idx] = sort(x);
		plot(x,y(idx),'.-','Color',c(j,:),'MarkerSize',20)
		y(x<0) = []; x(x<0) = [];
		% if length(x) > 4
		% 	ff = fit(x,y,ft,'StartPoint',[1 0 mean(x) 1],'Lower',[0 -10 0 0],'Upper',[5 5 max(x) 2]);
		% 	plot(logspace(-2,2,100),ff(logspace(-2,1,100)),'Color',c(j,:))
		% end
	end
	xlabel('S')
	ylabel('R')
	set(gca,'XScale','log','XLim',[1e-4 10],'XTick',[1e-4 1e-3 1e-2 1e-1 1 10])

	% plot the change in response vs. stimulus
	subplot(2,2,4); hold on
	for j = 1:length(all_back_levels)
		plot_these = sensillum == all_orns(i) & back_level == all_back_levels(j);
		x = dS(plot_these)+S_baseline(plot_these);
		y = dR(plot_these);
		[x,idx] = sort(x);
		plot(x,y(idx),'.-','Color',c(j,:),'MarkerSize',20)
		% y(x<0) = []; x(x<0) = [];
		% if length(x) > 4
		% 	ff = fit(x,y,ft,'StartPoint',[1 0 mean(x) 1],'Lower',[0 -10 0 0],'Upper',[5 5 max(x) 2]);
		% 	plot(logspace(-2,2,100),ff(logspace(-2,1,100)),'Color',c(j,:))
		% end

	end
	xlabel('S')
	ylabel('\DeltaR')
	set(gca,'XScale','log','XLim',[1e-4 10],'XTick',[1e-4 1e-3 1e-2 1e-1 1 10])


	prettyFig()

	labelFigure;

	if being_published	
		snapnow	
		delete(gcf)
	end
end


%% Kinetics 
% Now I look at the kinetics of the response. I fit exponentials to the rising and falling phases of the LFP to estimate timescales, and plot these timescales as functions of various things. (a) shows all the estimated timescales, plotted vs. the $r^2$ value of each estimate. I only retain estimates  with $r^2 > 0.8$ in subsequent plots. Colors as before, indicating background stimulus. Note that in (d), the on timescales increase steadily with background stimulus. 

% measure onset and offset rates for every trial
tau_on = NaN*sensillum;
tau_off = NaN*sensillum;
tau_on_r2 =  NaN*sensillum;
tau_off_r2 = NaN*sensillum;

ft_on = fittype('1 - exp(-x./tau)');
ft_off = fittype('exp(-x./tau)');
for i = 1:length(sensillum)
	if dR(i) < 0
		[~,loc] = min(R(1e3:1.5e3,i));
		x = R(1090:1e3+loc,i);
		x = -x; x = x-min(x); x = x/max(x);
		y = R(1e3+loc:end,i);
		y = -y; y = y - mean(y(end-5e3:end)); y = y/max(y);
	else
		[~,loc] = max(R(1e3:1.5e3,i));
		x = R(1090:1e3+loc,i);
		x = x-min(x); x = x/max(x);
		y = R(1e3+loc:end,i);
		y = y - mean(y(end-5e3:end)); y = y/max(y);
	end
	t = 1:length(x); t=t(:);
	try
		[ff,gof] = fit(t,x(:),ft_on,'Lower',0,'Start',100);
		tau_on(i) = ff.tau;
		tau_on_r2(i) = gof.rsquare;
	catch
	end


	t = 1:length(y); t=t(:);
	try	
		[ff, gof] = fit(t,y(:),ft_off,'Lower',0,'Start',100);
		tau_off(i) = ff.tau;
		tau_off_r2(i) = gof.rsquare;
	catch
	end
end

all_back_levels = unique(back_level);
c = parula(length(all_back_levels)+1);
figure('outerposition',[0 0 1311 902],'PaperUnits','points','PaperSize',[1311 902]); hold on
subplot(2,3,1); hold on
for j = 1:length(all_back_levels)
	plot_these = back_level == all_back_levels(j);
	plot(tau_on(plot_these),tau_on_r2(plot_these),'o','Color',c(j,:))
	plot(tau_off(plot_these),tau_off_r2(plot_these),'+','Color',c(j,:))
end
set(gca,'YLim',[0 1],'XScale','log','XLim',[10 1e4],'XTick',[10 100 1e3 1e4])
clear l
l(1) = plot(NaN,NaN,'ko');
l(2) = plot(NaN,NaN,'k+');
legend(l,{'\tau_{on}','\tau_{off}'},'Location','southeast')
xlabel('\tau (ms)')
ylabel('r^2')
plot([1e-3 1e5],[.8 .8],'k--')

subplot(2,3,2); hold on
for j = 1:length(all_back_levels)
	plot_these = back_level == all_back_levels(j) & tau_off_r2 > .8 & tau_on_r2 > .8;
	plot(tau_on(plot_these),tau_off(plot_these),'.','Color',c(j,:),'MarkerSize',20)
end
set(gca,'XLim',[0 500],'YLim',[0 1300])
xlabel('\tau_{on} (ms)')
ylabel('\tau_{off} (ms)')

subplot(2,3,3); hold on
for j = 1:length(all_back_levels)
	plot_these = back_level == all_back_levels(j) & tau_on_r2 > .8;
	plot(dR(plot_these),tau_on(plot_these),'o','Color',c(j,:))
	plot_these = back_level == all_back_levels(j) & tau_off_r2 > .8;
	plot(dR(plot_these),tau_off(plot_these),'+','Color',c(j,:))
end

l(1) = plot(NaN,NaN,'ko');
l(2) = plot(NaN,NaN,'k+');
legend(l,{'\tau_{on}','\tau_{off}'},'Location','northeast')
xlabel('\DeltaR')
ylabel('\tau (ms)')

subplot(2,3,4); hold on
for j = 1:length(all_back_levels)
	plot_these = back_level == all_back_levels(j) & tau_on_r2 > .8;
	plot(S_baseline(plot_these),tau_on(plot_these),'o','Color',c(j,:))
	plot_these = back_level == all_back_levels(j) & tau_off_r2 > .8;
	plot(S_baseline(plot_these),tau_off(plot_these),'+','Color',c(j,:))
end

l(1) = plot(NaN,NaN,'ko');
l(2) = plot(NaN,NaN,'k+');
legend(l,{'\tau_{on}','\tau_{off}'},'Location','northeast')
xlabel('S_{back}')
ylabel('\tau (ms)')


subplot(2,3,5); hold on
for j = 1:length(all_back_levels)
	plot_these = back_level == all_back_levels(j) & tau_on_r2 > .8;
	plot(R_baseline(plot_these),tau_on(plot_these),'o','Color',c(j,:))
	plot_these = back_level == all_back_levels(j) & tau_off_r2 > .8;
	plot(R_baseline(plot_these),tau_off(plot_these),'+','Color',c(j,:))
end

l(1) = plot(NaN,NaN,'ko');
l(2) = plot(NaN,NaN,'k+');
legend(l,{'\tau_{on}','\tau_{off}'},'Location','northeast')
xlabel('R_{baseline}')
ylabel('\tau (ms)')

subplot(2,3,6); hold on
for j = 1:length(all_back_levels)
	plot_these = back_level == all_back_levels(j) & tau_on_r2 > .8;
	plot(S_baseline(plot_these)+dS(plot_these),tau_on(plot_these),'o','Color',c(j,:))
	plot_these = back_level == all_back_levels(j) & tau_off_r2 > .8;
	plot(S_baseline(plot_these)+dS(plot_these),tau_off(plot_these),'+','Color',c(j,:))
end

l(1) = plot(NaN,NaN,'ko');
l(2) = plot(NaN,NaN,'k+');
legend(l,{'\tau_{on}','\tau_{off}'},'Location','northwest')
xlabel('S')
ylabel('\tau (ms)')

prettyFig();

labelFigure

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;


