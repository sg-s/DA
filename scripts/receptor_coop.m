
pHeader;


%% Receptor co-operativity 
% Are receptors co-operative? If we parameterise the input nonlinearity in a ORN using a Hill function, what is the $n$ parameter? 

%% 
% In the following figure, I plot three dose-response curves from two different odourant-receptor combinations, for both LFP and firing rate. (a) comes from a actual dose-response experiment that Carlotta did. (b-c) come from naturalistic odorant stimulation, where the stimulus was so broadly distributed that it spanned the entire sensitive range of the neuron.  

% get ab3A dose response
load(getPath(dataManager,'9e2a3e8664fe4f299b7ec435a0d2c2b1'));
v2struct(consolidated_data)
clear LFP Opt A_spikes* B_spikes*


figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on
x = []; y = [];
for i = 2:10
	this_resp = max(fA(1100:1500,paradigm==i));
	gas_dil = max(PID(1100:1500,paradigm==i));
	rm_this = gas_dil == 0 | this_resp == 0;
	gas_dil(rm_this) = []; this_resp(rm_this) = [];
	x = [x gas_dil];
	y = [y this_resp];
end

plot(x,y,'k+')
xlabel('ethyl acetate (V)')
ylabel('ab3A firing rate')
set(gca,'XScale','log')

ft = fittype(' hillFit(x,A,k,n,x_offset)');
ff = fit(x(:),y(:),ft,'StartPoint',[250 .1 1 0],'Lower',[1 1e-3 .1 0]);
l  =plot(sort(x),ff(sort(x)),'r');
legend(l,['n = ' oval(ff.n)],'Location','southeast')

% get the ab2A nat. stim data. 
cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);
whiff_stats = plotScaledNatStimWhiffStats(data,false);
x = vertcat(whiff_stats(2,:).stim_peaks);
y = vertcat(whiff_stats(2,:).peak_firing_rate);
rm_this = isnan(x) | isnan (y);
x(rm_this) = []; y(rm_this) = [];
subplot(1,3,2); hold on
plot(x,y,'k+')
ff = fit(x(:),y(:),ft,'StartPoint',[250 .1 1 0],'Lower',[1 1e-3 .1 0]);
l  =plot(sort(x),ff(sort(x)),'r');
legend(l,['n = ' oval(ff.n)],'Location','southeast')
set(gca,'XScale','log','XLim',[1e-2 10])
xlabel('2-butanone (V)')
ylabel('ab2A firing rate (Hz)')

x = vertcat(whiff_stats(2,:).stim_peaks);
y = -vertcat(whiff_stats(2,:).peak_LFP);
rm_this = isnan(x) | isnan (y);
x(rm_this) = []; y(rm_this) = [];
subplot(1,3,3); hold on
plot(x,y,'k+')
ff = fit(x(:),y(:),ft,'StartPoint',[25 .1 1 0],'Lower',[1 1e-3 .1 0]);
l  =plot(sort(x),ff(sort(x)),'r');
legend(l,['n = ' oval(ff.n)],'Location','southeast')
set(gca,'XScale','log','XLim',[1e-2 10])
xlabel('2-butanone (V)')
ylabel('ab2 LFP (mV)')


prettyFig()

labelFigure

if being_published	
	snapnow	
	delete(gcf)
end

%% Is this the real n?
% From the previous plot, it looks like $n<1$. I don't consider the kinetics of the response in the previous analysis. Is it possible that if we consider the kinetics, then we get a different n? How do we even determine this? One way is to fit an non-adapting kinetic model to the naturalistic stimulus data, and see what the best-fit $n$ is. The model is a simple binding-unbinding system, where the co-operativity of binding is allowed to vary. Despite initializing from $n=4$, best-fit models always had $n=1$. 

clear fd
for i = 1:3
	fd(i).stimulus = data(2).S(:,i);
	fd(i).response = data(2).X(:,i);
end

clear p
p.      n = 1;
p. k_plus = 196;
p.k_minus = 28;
p.      A = -18.7762;
p.      B = 0.5462;


%%
% In the following figure, I plot the LFP responses together with the best-fit model responses. Note that they do an OK job. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
time = 1e-3*(1:length(data(2).S));
clear l
l(1) = plot(time,data(2).X(:,1),'k');
[R, a] = KineticModelv1(data(2).S(:,1),p);
l(2) = plot(time,R,'r');
xlabel('Time (s)')
ylabel('LFP (mV)')
legend(l,{'ab2',['kinetic model, r^2 = ' oval(rsquare(R,data(2).X(:,1)))]},'Location','southwest')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% What if I force this model to have $n=4$, generate responses from this, and treat this as data, and then fit a Hill function to the responses to whiffs? 

figure('outerposition',[0 0 1502 500],'PaperUnits','points','PaperSize',[1502 500]); hold on
all_n = [2 3 4];

for ni = 1:length(all_n)
	n = all_n(ni);

	fake_data = data(2);
	p.n = n;
	for i = 1:3
		fake_data.X(:,i) = KineticModelv1(fake_data.S(:,i),p);
	end

	whiff_stats = plotScaledNatStimWhiffStats(fake_data,false);


	x = vertcat(whiff_stats(:).stim_peaks);
	y = -vertcat(whiff_stats(:).peak_LFP);
	rm_this = isnan(x) | isnan (y);
	x(rm_this) = []; y(rm_this) = [];
	subplot(1,3,ni); hold on
	plot(x,y,'k+')
	ff = fit(x(:),y(:),ft,'StartPoint',[25 .1 n 0],'Lower',[1 1e-3 .1 0]);
	l  =plot(sort(x),ff(sort(x)),'r');
	legend(l,['n = ' oval(ff.n)],'Location','southeast')
	set(gca,'XScale','log','XLim',[1e-2 10])
	xlabel('Stimulus (V)')
	ylabel('model LFP (mV)')
	title(['model n = ' oval(n)])
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Based on these simualtions, it looks like there is nothing wrong with my estimation of $n$ from the naturalistic stimulus data, and presumably from the dose response -- this is the "true" n. This leaves with our original problem: sparse stimuli suggest that $n = 1$, and the variance switching experiment suggests that $n=4$

%% Version Info
%
pFooter;


