
pHeader;


%% A closer look at lags
% In this document, we look at the effect of the mean stimulus on timescale in a little more detail. All the data comes from the Weber-Fechner experiment where we present fluctuating odor stimuli on top of backgrounds. 

a = 25e3;
z = 45e3;

clear cdata
cdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
cdata = cleanMSGdata(cdata);
cdata = rmfield(cdata,'K1');
cdata = rmfield(cdata,'K2');
v2struct(cdata)
time = 1e-3*(1:length(PID));

%% 1. Raw data: slow down in LFP
% In the following figure, we show the slowdown in the LFP in the raw traces by plotting the onset and offset responses of the LFP for various mean stimuli. Colors indicate mean stimuli, yellow = highest stimuli. 


% plot onset and offset responses of LFP and firing rate.
figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on
xlabel('Time (s)')
title('LFP responses to odor onset')
ylabel('LFP change')
subplot(1,3,2); hold on
title('LFP responses to odor offset')
xlabel('Time (s)')
ylabel('LFP change')

subplot(1,3,3); hold on
title('LFP responses to odor offset')
xlabel('Time (s)')
ylabel('LFP change (norm)')

c = parula(11);
sa = 4e3; sz = 7e3;
za = 54e3; zz = 57e3;
for i = 1:max(paradigm)
	X = LFP(:,paradigm==i);
	for j = 1:width(X)
		X(:,j) = X(:,j) - mean(X(1:5e3,j));
	end
	subplot(1,3,1); hold on
	plot(time(sa:sz),nanmean(X(sa:sz,:),2),'Color',c(i,:))

	X = LFP(:,paradigm==i);
	for j = 1:width(X)
		X(:,j) = X(:,j) - mean(X(50e3:55e3,j));
	end

	subplot(1,3,2); hold on
	plot(time(za:zz),nanmean(X(za:zz,:),2),'Color',c(i,:))

	subplot(1,3,3); hold on
	plot(time(za:zz),nanmean(X(za:zz,:),2)/max(nanmean(X(za:zz,:),2)),'Color',c(i,:))

end

prettyFig;


if being_published	
	snapnow	
	delete(gcf)
end

%% 2. Raw data: LFP kinetics vs. firing rate kinetics.
% Now, we compare LFP kinetics to firing rate kinetics to see if one is different from the other. IN the following figure, I'm plotting LFP and firing rate offset responses to odor. The LFP is plotted in blue (and inverted), and the firing rates are plotted in orange. Note that the firing rates always shut down quickly, while the LFPs don't. The color of the X-axes in each plot indicates the mean stimulus (more yellow= higher stimulus). 

figure('outerposition',[0 0 1300 791],'PaperUnits','points','PaperSize',[1300 791]); hold on
clear ax
max_paradigm = 8;
for i = 1:max_paradigm
	ax(i) = autoPlot(max_paradigm,i); hold on
end

sa = 4e3; sz = 7e3;
za = 54e3; zz = 57e3;
c = parula(max_paradigm+1);
for i = 1:max_paradigm
	X = LFP(:,paradigm==i);
	for j = 1:width(X)
		X(:,j) = X(:,j) - mean(X(50e3:55e3,j));
	end
	F = fA(:,paradigm==i);
	for j = 1:width(F)
		F(:,j) = F(:,j) - mean(F(50e3:55e3,j));
	end
	ax2 = plotyy(ax(i),time(za:zz),-nanmean(X(za:zz,:),2),time(za:zz),nanmean(F(za:zz,:),2));
	ax2(1).YLim = [-10 2];
	ax(i).YTick = [-10:2:2];
	ax2(2).YLim = [-30 30];
	ax2(2).YTick = [-30:10:30];
	ax2(1).XColor = c(i,:);
end

prettyFig;


if being_published	
	snapnow	
	delete(gcf)
end

%% 3. Cross correlation functions
% Now we compute the cross correlation functions between the stimulus, the LFP and the firing rate, and plot them, normalizing by peak to visualize the location of the peak (top row), and without normalising (bottom row): 


xSF = NaN(40e3+1,width(PID));
xSL = NaN(40e3+1,width(PID));
xLF = NaN(40e3+1,width(PID));

xa = 19.7e3;
xz = 20.3e3;

t = (1:length(xSF))*1e-3 - 20;

for i = 1:width(PID)
	% compute LFP lags
	s = PID(a:z,i)-mean(PID(a:z,i)); s = s/std(s);
	x = LFP(a:z,i)-mean(LFP(a:z,i)); x = x/std(x); x = -x;
	r = fA(a:z,i)-mean(fA(a:z,i)); r = r/std(r);

	temp_var = xcorr(x,s);  
	xSL(:,i) = temp_var/(z-a);

	temp_var = xcorr(r,s);  
	xSF(:,i) = temp_var/(z-a);

	temp_var = xcorr(r,x);  
	xLF(:,i) = temp_var/(z-a);
end

figure('outerposition',[0 0 1501 801],'PaperUnits','points','PaperSize',[1501 801]); hold on

c = parula(9);
for i = 1:length(paradigm)
	if paradigm(i) < 9
		subplot(2,3,1); hold on
		plot(t(xa:xz)*1e-3,xSL(xa:xz,i)/max(xSL(xa:xz,i)),'Color',c(paradigm(i),:))
		subplot(2,3,4); hold on
		plot(t(xa:xz)*1e-3,xSL(xa:xz,i),'Color',c(paradigm(i),:))
	end
end
subplot(2,3,1); hold on
xlabel('Lag (s)')
title('Stimulus -> LFP')
ylabel('xcorr (norm)')
set(gca,'YLim',[-1 1])

subplot(2,3,4); hold on
xlabel('Lag (s)')
title('Stimulus -> LFP')
ylabel('xcorr')
set(gca,'YLim',[-1 1])


for i = 1:length(paradigm)
	if paradigm(i) < 9
		subplot(2,3,2); hold on
		plot(t(xa:xz)*1e-3,xLF(xa:xz,i)/max(xLF(xa:xz,i)),'Color',c(paradigm(i),:))
		subplot(2,3,5); hold on
		plot(t(xa:xz)*1e-3,xLF(xa:xz,i),'Color',c(paradigm(i),:))
	end
end
subplot(2,3,2); hold on
title('LFP -> Firing')
xlabel('Lag (s)')
ylabel('xcorr (norm)')
set(gca,'YLim',[-1 1])

subplot(2,3,5); hold on
title('LFP -> Firing')
xlabel('Lag (s)')
ylabel('xcorr')
set(gca,'YLim',[-1 1])


c = parula(8);

for i = 1:length(paradigm)
	if paradigm(i) < 8
		subplot(2,3,3); hold on
		plot(t(xa:xz)*1e-3,xSF(xa:xz,i)/max(xSF(xa:xz,i)),'Color',c(paradigm(i),:))
		subplot(2,3,6); hold on
		plot(t(xa:xz)*1e-3,xSF(xa:xz,i),'Color',c(paradigm(i),:))
	end
end
subplot(2,3,3); hold on
title('Stimulus -> Firing')
xlabel('Lag (s)')
ylabel('xcorr (norm)')
set(gca,'YLim',[-1 1])

subplot(2,3,6); hold on
title('Stimulus -> Firing')
xlabel('Lag (s)')
ylabel('xcorr ')
set(gca,'YLim',[-1 1])

prettyFig;


if being_published	
	snapnow	
	delete(gcf)
end

%% 4. Filters.
% now, we back out filters from the stimulus, LFP and the firing rate. 
rf = 1;

if ~exist('K1','var')
	K1 = NaN(700,max(paradigm));
	K2 = K1;
	K3 = K1;
	for i = 1:max(paradigm)
		X = LFP(:,paradigm==i);
		R = fA(:,paradigm==i);
		S = PID(:,paradigm==i);

		temp_var = NaN(700,width(S));
		for j = 1:width(S)
			temp_var(:,j) = fitFilter2Data(nanmean(S(a:z,j),2),nanmean(R(a:z,j),2),'offset',100,'filter_length',700,'reg',rf);
		end
		K3(:,i) = nanmean(temp_var,2);

		temp_var = NaN(700,width(S));
		for j = 1:width(S)
			temp_var(:,j) = fitFilter2Data(nanmean(X(a:z,j),2),nanmean(R(a:z,j),2),'offset',100,'filter_length',700,'reg',rf);
		end
		K2(:,i) = nanmean(temp_var,2);

		temp_var = NaN(700,width(S));
		for j = 1:width(S)
			temp_var(:,j) = fitFilter2Data(nanmean(S(a:z,j),2),nanmean(X(a:z,j),2),'offset',100,'filter_length',700,'reg',rf);
		end
		K1(:,i) = nanmean(temp_var,2);
	end
end

figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,3); hold on
c = parula(9);
time = 1e-3*(1:length(fA));
filtertime = time(1:length(K3)) - .1;

for i = 1:max(paradigm)-2
	plot(filtertime,K3(:,i)/max(K3(50:200,i)),'Color',c(i,:))
end
set(gca,'YLim',[-1 1])
title('Stimulus -> Firing')
ylabel('Filter (norm)')

subplot(1,3,2); hold on
for i = 1:max(paradigm)-2
	plot(filtertime,K2(:,i)/min(K2(50:200,i)),'Color',c(i,:))
end
set(gca,'YLim',[-1 1])
title('LFP -> Firing')
ylabel('Filter (norm)')

subplot(1,3,1); hold on
for i = 1:max(paradigm)-2
	plot(filtertime,K1(:,i)/min(K1(50:200,i)),'Color',c(i,:))
end
set(gca,'YLim',[-1 1])
title('Stimulus -> LFP')
ylabel('Filter (norm)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% 
% Now, we predict LFP using the LFP filters, and then back out filters from the firing rate to the predicted LFP. 

% predict LFP
LFP_pred = NaN*LFP;
for i = 1:length(paradigm)
	LFP_pred(:,i) = convolve(time,PID(:,i),K1(:,paradigm(i)),filtertime);
end

% now back out filters from predicted LFP and firing rate

if ~exist('K4','var')
	K4 = NaN(700,max(paradigm));
	for i = 1:max(paradigm)
		X = LFP_pred(:,paradigm==i);
		R = fA(:,paradigm==i);

		temp_var = NaN(700,width(S));
		for j = 1:width(X)
			temp_var(:,j) = fitFilter2Data(nanmean(X(a:z,j),2),nanmean(R(a:z,j),2),'offset',100,'filter_length',700,'reg',rf);
		end
		K4(:,i) = nanmean(temp_var,2);
	end
end

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
c = parula(11);
for i = 1:max(paradigm)
	plot(filtertime,K4(:,i)/min(K4(50:200,i)),'Color',c(i,:))
end
set(gca,'YLim',[-1 1])
title('Predicted LFP -> Firing Rate')
ylabel('Filter (norm)')

prettyFig();


return


mean_stim = mean(PID(a:z,:));


% predictions using lowest LFP
f_syn = NaN(length(fA),10);
S = nanmean(LFP(:,paradigm==1),2);
time = 1e-3*(1:length(fA));
filtertime = time(1:length(K4)) - .1;
for i = 1:max(paradigm)
	f_syn(:,i) = convolve(time,S,K4(:,i),filtertime);
	% normalise
	f_syn(:,i) = f_syn(:,i) - mean(f_syn(a:z,i));
	f_syn(:,i) = f_syn(:,i)/std(f_syn(a:z,i));
end

figure, hold on
c = parula(11);
for i = 1:max(paradigm)
	plot(time,f_syn(:,i),'Color',c(i,:))
end





load(getPath(dataManager,'aeb361c027b71938021c12a6a12a85cd'),'-mat');


odour_thresh = .05;

example_orn = 6;

ab3.PID = [od(example_orn).stimulus];
for i = 1:width(ab3.PID)
	ab3.PID(:,i) = ab3.PID(:,i) -  mean(ab3.PID(1:5e3,i));
end


figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
% show the whiff durations 
ax(3) = subplot(1,2,1); hold on

whiff_durations = []; 
for i = 1:width(ab3.PID)
	[ons,offs] = computeOnsOffs(ab3.PID(:,i) > odour_thresh);
	whiff_durations =  [whiff_durations; offs-ons];
end
whiff_durations = nonzeros(whiff_durations);
[y,x] = histcounts(whiff_durations,50); x(1) = [];
y = y/sum(y);
a = 1; m = fittype('a + n*x');
xx = vectorise(log(x)); yy = vectorise(log(y));
ff = fit(xx(yy>-Inf),yy(yy>-Inf),m,'Upper',[Inf -1.5],'Lower',[-Inf -1.5],'StartPoint',[300 -1.5]);
plot(ax(3),x,y,'k+')
plot(ax(3),x,exp(ff(log(x))),'r')
ylabel(ax(3),'Probability')
set(ax(3),'YScale','log','XScale','log','XTick',[1e1 1e2 1e3 1e4],'XLim',[10 1.1e4])
xlabel(ax(3),'Whiff duration (ms)')


% show the blank durations 
ax(3) = subplot(1,2,2); hold on

whiff_durations = []; 
for i = 1:width(ab3.PID)
	[ons,offs] = computeOnsOffs(ab3.PID(:,i) < odour_thresh);
	whiff_durations =  [whiff_durations; offs-ons];
end
whiff_durations = nonzeros(whiff_durations);
[y,x] = histcounts(whiff_durations,50); x(1) = [];
y = y/sum(y);
a = 1; m = fittype('a + n*x');
xx = vectorise(log(x)); yy = vectorise(log(y));
ff = fit(xx(yy>-Inf),yy(yy>-Inf),m,'Upper',[Inf -1.5],'Lower',[-Inf -1.5],'StartPoint',[300 -1.5]);
plot(ax(3),x,y,'k+')
plot(ax(3),x,exp(ff(log(x))),'r')
ylabel(ax(3),'Probability')
set(ax(3),'YScale','log','XScale','log','XTick',[1e1 1e2 1e3 1e4],'XLim',[10 1.1e4])
xlabel(ax(3),'Blank duration (ms)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
%
pFooter;






