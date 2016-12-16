


pHeader;

%% Paired Pulses 
% In this document, I attempt to estimate the timescale of gain control using paired pulses, as in Cao et al. PNAS. 

%% ab3 LFP
% First, I plot the data we got from LFP responses in the ab3 sensillum to ethyl acetate odourant. In the following figure, I plot the stimulus and the LFP responses from 4 differnet ORNs, colour coded by the inter-pulse interval, similar to the plot in Cao et al. PNAS. 

   ;;;    ;;;;;;;;   ;;;;;;;     ;;       ;;;;;;;; ;;;;;;;;  
  ;; ;;   ;;     ;; ;;     ;;    ;;       ;;       ;;     ;; 
 ;;   ;;  ;;     ;;        ;;    ;;       ;;       ;;     ;; 
;;     ;; ;;;;;;;;   ;;;;;;;     ;;       ;;;;;;   ;;;;;;;;  
;;;;;;;;; ;;     ;;        ;;    ;;       ;;       ;;        
;;     ;; ;;     ;; ;;     ;;    ;;       ;;       ;;        
;;     ;; ;;;;;;;;   ;;;;;;;     ;;;;;;;; ;;       ;;        


% get the raw data
cdata  = consolidateDataPairedPulses('/Volumes/sgs_data/Dropbox (emonetlab)/users/srinivas_gs/data/DA-paper/data-for-paper/paired-pulses/ok/');
v2struct(cdata)

% remove the baseline from all LFP traces
for i = 1:length(orn)
	LFP(:,i) = LFP(:,i) - mean(LFP(4e3:5e3,i));
end

% remove the baseline from all PID traces
for i = 1:length(orn)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

time = 1e-3*(1:length(PID));



c = parula(length(unique(pulse_seperation))+1);
all_pulse_sep = unique(pulse_seperation);


figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,1,1); hold on
for i = length(all_pulse_sep):-1:1
	plot(time,mean(PID(:,pulse_seperation == all_pulse_sep(i)),2),'Color',c(i,:));
end
xlabel('Time (s)')
ylabel('Stimulus (V)')
set(gca,'XLim',[0 20])

subplot(2,1,2); hold on
for i = length(all_pulse_sep):-1:1
	plot(time,mean(LFP(:,pulse_seperation == all_pulse_sep(i)),2),'Color',c(i,:));
end
xlabel('Time (s)')
ylabel('\DeltaLFP (mV)')
set(gca,'XLim',[0 20])


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% I now plot the ratio of the response to the second pulse to the first as a function of the pulse spacing. There are many various to do this, and I do all of them. First, in (a), I plot the ratios of peak response to the second pulse to the response to the first pulse, as in Cao et al. Our data differs from theirs: we never see a decrease in the peak response amplitude. It instead increases at short pulse separations, before decaying to 1. Next, in (b), I plot the ratio of the increase in response to the two pulses. Finally, in (c), I plot the ratio of the gains measured at the two pulses, measured as change in response divided by change in stimulus. In (b) and (c), I fit an exponential to the data to determine the timescale of putative gain control. Note that our timescales are much faster than what was reported in Cao et al., and agree with the timescales we reported in the paper (<1s). 


% measure peak response for each pulse
B1 = NaN*orn;
B2 = NaN*orn;
P1 = NaN*orn;
P2 = NaN*orn;
S1 = NaN*orn;
S2 = NaN*orn;

for i = 1:length(orn)
	x = LFP(:,i);
	B1(i) = mean(x(4e3:5e3));

	s = PID(:,i);


	% first pulse peak
	a = 5.05e3; z = a + 150;
	[~,loc_A] = min(x(a:z));
	loc_A = round(a + loc_A);
	P1(i) = x(loc_A);

	[S1(i),sA] = max(s(a:z));
	sA = sA + a;

	% second pulse peak
	a = round(5.05e3 + pulse_seperation(i)*1e3 + 50); z = a + 150;
	[~,loc_B] = min(x(a:z));
	loc_B = round(a + loc_B);
	P2(i) = x(loc_B);

	[S2(i),sB] = max(s(a:z));
	sB = sB + a;

	% find the maximum value b/w the peaks
	[~,loc_C] = max(x(loc_A:loc_B));
	loc_C = round(loc_C + loc_A);

	B2(i) = x(loc_C);

	% figure, hold on
	% plot(s,'k')
	% plot(sA,s(sA),'ro')
	% plot(sB,s(sB),'bo')

	% plot(loc_B,x(loc_B),'bo')
	% plot(loc_C,x(loc_C),'go')
	% plot(loc_A,x(loc_A),'ro')
end


figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on

subplot(1,3,1); hold on
plot(pulse_seperation,-P2./-P1,'k+')
xlabel('\Delta T (s)')
ylabel('R_{2}/R_0')

subplot(1,3,2); hold on
y = -(P2-B2)./-(P1-B1);
plot(pulse_seperation,y,'k+')
xlabel('\Delta T (s)')
ylabel('\DeltaR_2/\DeltaR_1')


cf = fittype('1-exp(-x./tau)');
ff = fit(pulse_seperation,y(:),cf,'Start',.5,'Lower',0,'Upper',10);
x = logspace(-1,1,1e3);
l = plot(x,ff(x),'r');

legend(l,['\tau = ' oval(ff.tau), 's r^2 = ' oval(rsquare(y,ff(pulse_seperation)))],'Location','southeast')


subplot(1,3,3); hold on
y2 = -(P2-B2)./(S2);
y1 = -(P1-B1)./(S1);
y = y2./y1;
plot(pulse_seperation,y,'k+')
xlabel('\Delta T (s)')
ylabel('G_2/G_1')

cf = fittype('1-exp(-x./tau)');
ff = fit(pulse_seperation,y(:),cf,'Start',.5,'Lower',0,'Upper',10);
x = logspace(-1,1,1e3);
l = plot(x,ff(x),'r');

legend(l,['\tau = ' oval(ff.tau), 's r^2 = ' oval(rsquare(y,ff(pulse_seperation)))],'Location','southeast')

prettyFig();

labelFigure()

if being_published
	snapnow
	delete(gcf)
end


   ;;;    ;;;;;;;;   ;;;;;;;     ;;;;;;;; ;;;; ;;;;;;;;  ;;;; ;;    ;;  ;;;;;;   
  ;; ;;   ;;     ;; ;;     ;;    ;;        ;;  ;;     ;;  ;;  ;;;   ;; ;;    ;;  
 ;;   ;;  ;;     ;;        ;;    ;;        ;;  ;;     ;;  ;;  ;;;;  ;; ;;        
;;     ;; ;;;;;;;;   ;;;;;;;     ;;;;;;    ;;  ;;;;;;;;   ;;  ;; ;; ;; ;;   ;;;; 
;;;;;;;;; ;;     ;;        ;;    ;;        ;;  ;;   ;;    ;;  ;;  ;;;; ;;    ;;  
;;     ;; ;;     ;; ;;     ;;    ;;        ;;  ;;    ;;   ;;  ;;   ;;; ;;    ;;  
;;     ;; ;;;;;;;;   ;;;;;;;     ;;       ;;;; ;;     ;; ;;;; ;;    ;;  ;;;;;;   


%% ab3A firing rate responses
% Now I make the same plots for firing rate of the ab3A ORN. 


c = parula(length(unique(pulse_seperation))+1);
all_pulse_sep = unique(pulse_seperation);

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,1,1); hold on
for i = length(all_pulse_sep):-1:1
	plot(time,mean(PID(:,pulse_seperation == all_pulse_sep(i)),2),'Color',c(i,:));
end
xlabel('Time (s)')
ylabel('Stimulus (V)')
set(gca,'XLim',[0 20])

subplot(2,1,2); hold on
for i = length(all_pulse_sep):-1:1
	temp = fA(:,pulse_seperation == all_pulse_sep(i));
	temp(:,sum(temp)==0) = [];
	plot(time,mean(temp,2),'Color',c(i,:));
end
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
set(gca,'XLim',[0 20])


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% I now plot the ratio of the response to the second pulse to the first as a function of the pulse spacing. There are many various to do this, and I do all of them. First, in (a), I plot the ratios of peak response to the second pulse to the response to the first pulse, as in Cao et al. Our data differs from theirs: we never see a decrease in the peak response amplitude. It instead increases at short pulse separations, before decaying to 1. Next, in (b), I plot the ratio of the increase in response to the two pulses. Finally, in (c), I plot the ratio of the gains measured at the two pulses, measured as change in response divided by change in stimulus. In (b) and (c), I fit an exponential to the data to determine the timescale of putative gain control. Note that our timescales are much faster than what was reported in Cao et al., and agree with the timescales we reported in the paper (<1s). 

% measure peak response for each pulse
B1 = NaN*orn;
B2 = NaN*orn;
P1 = NaN*orn;
P2 = NaN*orn;
S1 = NaN*orn;
S2 = NaN*orn;

for i = 1:length(orn)
	x = fA(:,i);
	B1(i) = mean(x(4e3:5e3));

	s = PID(:,i);


	% first pulse peak
	a = 5.05e3; z = a + 150;
	[~,loc_A] = max(x(a:z));
	loc_A = round(a + loc_A);
	P1(i) = x(loc_A);

	[S1(i),sA] = max(s(a:z));
	sA = sA + a;

	% second pulse peak
	a = round(5.05e3 + pulse_seperation(i)*1e3 + 50); z = a + 150;
	[~,loc_B] = max(x(a:z));
	loc_B = round(a + loc_B);
	P2(i) = x(loc_B);

	[S2(i),sB] = max(s(a:z));
	sB = sB + a;

	% find the min value b/w the peaks
	[~,loc_C] = min(x(loc_A:loc_B));
	loc_C = round(loc_C + loc_A);

	B2(i) = x(loc_C);

end

rm_this = P1 == 0 | P2 == 0;
P1(rm_this) = NaN;
P2(rm_this) = NaN;
B1(rm_this) = NaN;
B2(rm_this) = NaN;
S1(rm_this) = NaN;
S2(rm_this) = NaN;

figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on

subplot(1,3,1); hold on
plot(pulse_seperation,-P2./-P1,'k+')
xlabel('\Delta T (s)')
ylabel('R_{a}/R_0')

subplot(1,3,2); hold on
y = -(P2-B2)./-(P1-B1);
plot(pulse_seperation,y,'k+')
xlabel('\Delta T (s)')
ylabel('\DeltaR_{2}/\DeltaR_1')


cf = fittype('1-exp(-x./tau)');

ff = fit(pulse_seperation(~isnan(y)),nonnans(y(:)),cf,'Start',.5,'Lower',0,'Upper',10);
x = logspace(-1,1,1e3);
l = plot(x,ff(x),'r');

legend(l,['\tau = ' oval(ff.tau), 's r^2 = ' oval(rsquare(y,ff(pulse_seperation)))],'Location','southeast')


subplot(1,3,3); hold on
y2 = -(P2-B2)./(S2);
y1 = -(P1-B1)./(S1);
y = y2./y1;
plot(pulse_seperation,y,'k+')
xlabel('\Delta T (s)')
ylabel('G_2/G_1')

cf = fittype('1-exp(-x./tau)');
ff = fit(pulse_seperation(~isnan(y)),nonnans(y(:)),cf,'Start',.5,'Lower',0,'Upper',10);
x = logspace(-1,1,1e3);
l = plot(x,ff(x),'r');

legend(l,['\tau = ' oval(ff.tau), 's r^2 = ' oval(rsquare(y,ff(pulse_seperation)))],'Location','southeast')

prettyFig();

labelFigure

if being_published
	snapnow
	delete(gcf)
end

   ;;;    ;;;;;;;;   ;;;;;;;     ;;       ;;;;;;;; ;;;;;;;;  
  ;; ;;   ;;     ;; ;;     ;;    ;;       ;;       ;;     ;; 
 ;;   ;;  ;;     ;;        ;;    ;;       ;;       ;;     ;; 
;;     ;; ;;;;;;;;   ;;;;;;;     ;;       ;;;;;;   ;;;;;;;;  
;;;;;;;;; ;;     ;; ;;           ;;       ;;       ;;        
;;     ;; ;;     ;; ;;           ;;       ;;       ;;        
;;     ;; ;;;;;;;;  ;;;;;;;;;    ;;;;;;;; ;;       ;;        


%% ab2 LFP responses
% That doesn't look very convincing either. What if we do the same experiment on ab2? The idea is that since ab2 is much more sensitive to this odour, these pulses will be saturating to the neuron, and thus more likely to induce adaptation. 

cdata  = consolidateDataPairedPulses('/Volumes/sgs_data/Dropbox (emonetlab)/users/srinivas_gs/data/DA-paper/data-for-paper/paired-pulses/ab2/');
v2struct(cdata);

% remove the baseline from all LFP traces
for i = 1:length(orn)
	LFP(:,i) = LFP(:,i) - mean(LFP(4e3:5e3,i));
end

% remove the baseline from all PID traces
for i = 1:length(orn)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

time = 1e-3*(1:length(PID));


c = parula(length(unique(pulse_seperation))+1);
all_pulse_sep = unique(pulse_seperation);


figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,1,1); hold on
for i = length(all_pulse_sep):-1:1
	plot(time,mean(PID(:,pulse_seperation == all_pulse_sep(i)),2),'Color',c(i,:));
end
xlabel('Time (s)')
ylabel('Stimulus (V)')
set(gca,'XLim',[0 20])

subplot(2,1,2); hold on
for i = length(all_pulse_sep):-1:1
	plot(time,mean(LFP(:,pulse_seperation == all_pulse_sep(i)),2),'Color',c(i,:));
end
xlabel('Time (s)')
ylabel('ab2 \DeltaLFP (mV)')
set(gca,'XLim',[0 20])


prettyFig()

labelFigure

if being_published	
	snapnow	
	delete(gcf)
end



%%
% I now plot the ratio of the response to the second pulse to the first as a function of the pulse spacing. There are many various to do this, and I do all of them. 


% measure peak response for each pulse
B1 = NaN*orn;
B2 = NaN*orn;
P1 = NaN*orn;
P2 = NaN*orn;
S1 = NaN*orn;
S2 = NaN*orn;

for i = 1:length(orn)
	x = LFP(:,i);
	B1(i) = mean(x(4e3:5e3));

	s = PID(:,i);


	% first pulse peak
	a = 5.05e3; z = a + 150;
	[~,loc_A] = min(x(a:z));
	loc_A = round(a + loc_A);
	P1(i) = x(loc_A);

	[S1(i),sA] = max(s(a:z));
	sA = sA + a;

	% second pulse peak
	a = round(5.05e3 + pulse_seperation(i)*1e3 + 50); z = a + 150;
	[~,loc_B] = min(x(a:z));
	loc_B = round(a + loc_B);
	P2(i) = x(loc_B);

	[S2(i),sB] = max(s(a:z));
	sB = sB + a;

	% find the maximum value b/w the peaks
	[~,loc_C] = max(x(loc_A:loc_B));
	loc_C = round(loc_C + loc_A);

	B2(i) = x(loc_C);

	% figure, hold on
	% plot(s,'k')
	% plot(sA,s(sA),'ro')
	% plot(sB,s(sB),'bo')

	% plot(loc_B,x(loc_B),'bo')
	% plot(loc_C,x(loc_C),'go')
	% plot(loc_A,x(loc_A),'ro')
end


figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on

subplot(1,3,1); hold on
plot(pulse_seperation,-P2./-P1,'k+')
xlabel('\Delta T (s)')
ylabel('R_{2}/R_0')

subplot(1,3,2); hold on
y = -(P2-B2)./-(P1-B1);
plot(pulse_seperation,y,'k+')
xlabel('\Delta T (s)')
ylabel('\DeltaR_2/\DeltaR_1')


cf = fittype('1-exp(-x./tau)');
ff = fit(pulse_seperation,y(:),cf,'Start',.5,'Lower',0,'Upper',10);
x = logspace(-1,1,1e3);
l = plot(x,ff(x),'r');

legend(l,['\tau = ' oval(ff.tau), 's r^2 = ' oval(rsquare(y,ff(pulse_seperation)))],'Location','southeast')


subplot(1,3,3); hold on
y2 = -(P2-B2)./(S2);
y1 = -(P1-B1)./(S1);
y = y2./y1;
plot(pulse_seperation,y,'k+')
xlabel('\Delta T (s)')
ylabel('G_2/G_1')

cf = fittype('1-exp(-x./tau)');
ff = fit(pulse_seperation,y(:),cf,'Start',.5,'Lower',0,'Upper',10);
x = logspace(-1,1,1e3);
l = plot(x,ff(x),'r');

legend(l,['\tau = ' oval(ff.tau), 's r^2 = ' oval(rsquare(y,ff(pulse_seperation)))],'Location','southeast')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end





%% Version Info
%
pFooter;


