


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
cdata  = consolidateDataPairedPulses('/data/DA-paper/data-for-paper/paired-pulses/ok/');
v2struct(cdata)

% remove the baseline from all LFP traces
for i = 1:length(orn)
	LFP(:,i) = LFP(:,i) - mean(LFP(4e3:5e3,i));
end

% remove the baseline from all PID traces
for i = 1:length(orn)
	PID(:,i) = PID(:,i) - min(PID(1:5e3,i));
end

time = 1e-3*(1:length(PID));



c = parula(length(unique(pulse_seperation))+1);
all_pulse_sep = unique(pulse_seperation);


figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(2,1,1); hold on
clear data
for i = length(all_pulse_sep):-1:1
	data(i).stimulus = mean(PID(:,pulse_seperation == all_pulse_sep(i)),2);
	plot(time,data(i).stimulus,'Color',c(i,:));
end
xlabel('Time (s)')
ylabel('Stimulus (V)')
set(gca,'XLim',[0 20])

subplot(2,1,2); hold on
for i = length(all_pulse_sep):-1:1
	data(i).response = mean(LFP(:,pulse_seperation == all_pulse_sep(i)),2);
	plot(time,data(i).response,'Color',c(i,:));
	data(i).response = -data(i).response;
end
xlabel('Time (s)')
ylabel('\DeltaLFP (mV)')
set(gca,'XLim',[0 20])

suptitle('ab3 LFP')
prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% I now plot the ratio of the response to the second pulse to the first as a function of the pulse spacing. There are many various to do this, and I do all of them. First, in (a), I plot the ratios of peak response to the second pulse to the response to the first pulse, as in Cao et al. Our data differs from theirs: we never see a decrease in the peak response amplitude. It instead increases at short pulse separations, before decaying to 1.  Finally, in (b), I plot the ratio of the gains measured at the two pulses, measured as change in response divided by change in stimulus. In (b), I fit an exponential to the data to determine the timescale of putative gain control. Note that our timescales are much faster than what was reported in Cao et al., and agree with the timescales we reported in the paper (<1s). 


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
ylabel('R_{2}/R_1')

subplot(1,3,2); hold on
plot(pulse_seperation,S2./S1,'k+')
xlabel('\Delta T (s)')
ylabel('S_2/S_1')


% cf = fittype('1-exp(-x./tau)');
% ff = fit(pulse_seperation,y(:),cf,'Start',.5,'Lower',0,'Upper',10);
% x = logspace(-1,1,1e3);
% l = plot(x,ff(x),'r');

% legend(l,['\tau = ' oval(ff.tau), 's r^2 = ' oval(rsquare(y,ff(pulse_seperation)))],'Location','southeast')


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

suptitle('ab3 LFP')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

 ;;;;;;  ;;;; ;;     ;; ;;     ;; ;;          ;;;    ;;;;;;;; ;;;;  ;;;;;;;  ;;    ;;  ;;;;;;  
;;    ;;  ;;  ;;;   ;;; ;;     ;; ;;         ;; ;;      ;;     ;;  ;;     ;; ;;;   ;; ;;    ;; 
;;        ;;  ;;;; ;;;; ;;     ;; ;;        ;;   ;;     ;;     ;;  ;;     ;; ;;;;  ;; ;;       
 ;;;;;;   ;;  ;; ;;; ;; ;;     ;; ;;       ;;     ;;    ;;     ;;  ;;     ;; ;; ;; ;;  ;;;;;;  
      ;;  ;;  ;;     ;; ;;     ;; ;;       ;;;;;;;;;    ;;     ;;  ;;     ;; ;;  ;;;;       ;; 
;;    ;;  ;;  ;;     ;; ;;     ;; ;;       ;;     ;;    ;;     ;;  ;;     ;; ;;   ;;; ;;    ;; 
 ;;;;;;  ;;;; ;;     ;;  ;;;;;;;  ;;;;;;;; ;;     ;;    ;;    ;;;;  ;;;;;;;  ;;    ;;  ;;;;;;  


%% Simulations with NL model
% To check whether any of this makes sense, I did some simulations with a simple NL model. To choose non-pathalogical parameters, I fit the NL model to this data (the LFP responses). and then used that NL model to generate synthetic data using the real stimulus. I treated this as data, and repeated the analysis as before. From the plots below, it looks like the ratio of gains (right) does not relax towards 1 as an exponential (cf. the data before). This is good -- the NL model should not show any sign of adaptation here. 

clear p
p.Hill_n = 1;
p.Hill_K = 1.1904;
p.  tau1 = 56.8594;
p.  tau2 = 400;
p.     n = 1;
p.     A = 0;
p.     C = 10.3040;

% generate responses using pNL
R = NaN*fA;
for i = 1:size(fA,2)
	R(:,i) = pNL(PID(:,i),p);
end
R = -R;

% measure peak response for each pulse
B1 = NaN*orn;
B2 = NaN*orn;
P1 = NaN*orn;
P2 = NaN*orn;
S1 = NaN*orn;
S2 = NaN*orn;

for i = 1:length(orn)
	x = R(:,i);
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

end


figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on

subplot(1,3,1); hold on
plot(pulse_seperation,-P2./-P1,'k+')
xlabel('\Delta T (s)')
ylabel('R_{2}/R_1')

subplot(1,3,2); hold on
plot(pulse_seperation,S2./S1,'k+')
xlabel('\Delta T (s)')
ylabel('S_2/S_1')


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

suptitle('NL model simulations')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Can this paired pulse protocol actually recover the timescale of gain control? To test this, I fit an adapting NL model to the data, where an adapting mechanism changes the $k_D$ of the input nonlinearity. I then set the timescale of this update to 500ms by hand, and see if I can recover this timescale from this analysis. As we see from the plots below, it fails spectacularly, and the recovered timescale is not even close to the actual timescale of adaptation. 

clear p
p.   k0 = 1.0337;
p.tau_z = 16.0938;
p.tau_z = 250;
p.    B = 9.3125;
p.  n_z = 2;
p.    n = 1.0781;
p. tau1 = 62.9375;
p. tau2 = 63.5312;
p.  n_y = 1;
p.    A = 0;
p.    C = 17.4375;

% generate responses using pNL
R = NaN*fA;
for i = 1:size(fA,2)
	R(:,i) = aNL(abs(PID(:,i)),p);
end
R = -R;

% measure peak response for each pulse
B1 = NaN*orn;
B2 = NaN*orn;
P1 = NaN*orn;
P2 = NaN*orn;
S1 = NaN*orn;
S2 = NaN*orn;

for i = 1:length(orn)
	x = R(:,i);
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

end


figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on

subplot(1,3,1); hold on
plot(pulse_seperation,-P2./-P1,'k+')
xlabel('\Delta T (s)')
ylabel('R_{2}/R_1')

subplot(1,3,2); hold on
plot(pulse_seperation,S2./S1,'k+')
xlabel('\Delta T (s)')
ylabel('S_2/S_1')


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

suptitle('adapting NL model simulations, \tau = 500ms')
prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% 
% So this is really not working. To figure out what's going on, I switch to an idealized stimulus that is a perfect square pulse (top panel). I then use the same model as before to generate responses to this stimulus (middle panel). We see the characteristic decrease in the peak response for short inter-pulse intervals, as in Cao et al. However, when I plot the ratio of the 2nd to 1st responses, I see something interesting: the ratio typically always goes above 1 for short timescales. It looks exponential only in the terminal segment. I fit that segment with an exponential, and extract a timescale as in Cao et al. I do this for three differnet simulations, with three different timescales of gain control (bottom three plots). Cao et al.'s method fails to recover the actual timescale of gain control in all three cases. For comparison, I plot a blue line at the location of the minimum of the ratio (a simple heuristic that intuitively is the timescale of gain control) -- and this does much better at estimating the true timescale of gain control than Cao's exponential fit. 

%%
% In conclusion, I show with a very simple adapting model that the paired-pulse protocol used by Cao and the analysis method that goes with it can be terribly misled. There is no guarantee that the timescale of the exponential fit has anything to do with the timescale of gain control. 

pulse_width = 50;
pulse_height = .6;
aps = round(1e3*unique(pulse_seperation));
S = zeros(length(PID),length(aps));
for i = 1:length(aps)
	S(5e3:5e3+pulse_width,i) = pulse_height;
	S(5e3+pulse_width+aps(i):5e3+2*pulse_width+aps(i),i) = pulse_height;
end

% generate responses
R = 0*S;
for i = 1:length(aps)
	 [R(:,i),Ky,Kz,k_D,x1] = aNL(S(:,i),p);
end


figure('outerposition',[0 0 901 888],'PaperUnits','points','PaperSize',[901 888]); hold on
subplot(3,1,1); hold on
set(gca,'ColorOrder',parula(11))
plot(S)
ylabel('Stimulus')

subplot(3,1,2); hold on
set(gca,'ColorOrder',parula(11))
plot(R)
ylabel('Response')
suptitle('adapting NL model simulations')
title('\tau_{adaptation} = 500ms')

aps = round(logspace(1,4,100));
S = zeros(length(PID),length(aps));
for i = 1:length(aps)
	S(5e3:5e3+pulse_width,i) = pulse_height;
	S(5e3+pulse_width+aps(i):5e3+2*pulse_width+aps(i),i) = pulse_height;
end

% do the analysis for three different timescales 
all_tau = [300 600 2000];
for k = 1:length(all_tau)
	p.tau_z = all_tau(k)/p.n_z;
	% generate responses
	R = 0*S;
	for i = 1:length(aps)
		R(:,i) = aNL(S(:,i),p);
	end

	% find peaks
	P1 = NaN*aps;
	P2 = NaN*aps;
	for i = 1:length(aps)
		[ons,offs] = computeOnsOffs(S(:,i));
		P1(i) = max(R(ons(1):ons(2),i));
		P2(i) = max(R(ons(2):end,i));
	end

	subplot(3,3,6+k); hold on
	plot(aps,P2./P1,'k.')
	x = aps(:);
	y = vectorise(P2./P1);
	cf = fittype('1-exp(-x./tau)');
	[~,loc]=min(P2./P1);
	ff = fit(x(loc:end),y(loc:end),cf,'Start',400,'Lower',0,'Upper',Inf);
	clear l
	l(1) = plot(aps,ff(aps),'r');
	title(['Actual \tau = ' oval(all_tau(k)) 'ms'])
	set(gca,'XScale','log','YLim',[0 4],'XTick',[10 1e2 1e3 1e4])
	l(2) = plot([aps(loc) aps(loc)],[0 2.5],'b');
	legend(l,{['\tau = ' oval(ff.tau) 'ms'],['\tau = ' oval(aps(loc)) 'ms']})
	xlabel('Pulse spacing (ms)')
	ylabel('R_2/R_1')
end

prettyFig('fs',15);

if being_published
	snapnow
	delete(gcf)
end


%% Addendum
% It now strikes me that this analysis method confuses the timescale of gain control with the amount of gain control. I did some simulations where I increased the strength of adaptation, but kept all other parameters the same. Peaks of subsequent responses now were smaller than the first one -- but Cao's method still failed to recover the actual gain control timescale. 


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

suptitle('ab3A firing rate')
prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% I now plot the ratio of the response to the second pulse to the first as a function of the pulse spacing. There are many various to do this, and I do all of them. First, in (a), I plot the ratios of peak response to the second pulse to the response to the first pulse, as in Cao et al. Our data differs from theirs: we never see a decrease in the peak response amplitude. It instead increases at short pulse separations, before decaying to 1.  Finally, in (b), I plot the ratio of the gains measured at the two pulses, measured as change in response divided by change in stimulus. In (b), I fit an exponential to the data to determine the timescale of putative gain control. Note that our timescales are much faster than what was reported in Cao et al., and agree with the timescales we reported in the paper (<1s). 

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

figure('outerposition',[0 0 1001 500],'PaperUnits','points','PaperSize',[1001 500]); hold on

subplot(1,2,1); hold on
plot(pulse_seperation,-P2./-P1,'k+')
xlabel('\Delta T (s)')
ylabel('R_2/R_1')

% subplot(1,3,2); hold on
% y = -(P2-B2)./-(P1-B1);
% plot(pulse_seperation,y,'k+')
% xlabel('\Delta T (s)')
% ylabel('\DeltaR_{2}/\DeltaR_1')


% cf = fittype('1-exp(-x./tau)');

% ff = fit(pulse_seperation(~isnan(y)),nonnans(y(:)),cf,'Start',.5,'Lower',0,'Upper',10);
% x = logspace(-1,1,1e3);
% l = plot(x,ff(x),'r');

% legend(l,['\tau = ' oval(ff.tau), 's r^2 = ' oval(rsquare(y,ff(pulse_seperation)))],'Location','southeast')


subplot(1,2,2); hold on
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

suptitle('ab3A firing rate')
prettyFig();

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

suptitle('ab2 LFP')
prettyFig()

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
ylabel('R_2/R_1')

subplot(1,3,2); hold on
plot(pulse_seperation,S2./S1,'k+')
xlabel('\Delta T (s)')
ylabel('S_2/S_1')


% cf = fittype('1-exp(-x./tau)');
% ff = fit(pulse_seperation,y(:),cf,'Start',.5,'Lower',0,'Upper',10);
% x = logspace(-1,1,1e3);
% l = plot(x,ff(x),'r');

% legend(l,['\tau = ' oval(ff.tau), 's r^2 = ' oval(rsquare(y,ff(pulse_seperation)))],'Location','southeast')


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

suptitle('ab2 LFP')
prettyFig();



if being_published
	snapnow
	delete(gcf)
end





%% Version Info
%
pFooter;


