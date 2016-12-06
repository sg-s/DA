
pHeader;

%% Can NL models typically lead to variance gain control?
% In this document, I check if a simple NL model can lead to something that we think of as variance gain control (along the lines of Nemenman's paper). 

;;      ;; ;;     ;; ;;;; ;;;;;;;; ;;;;;;;;    ;;    ;;  ;;;;;;;  ;;;;  ;;;;;;  ;;;;;;;; 
;;  ;;  ;; ;;     ;;  ;;     ;;    ;;          ;;;   ;; ;;     ;;  ;;  ;;    ;; ;;       
;;  ;;  ;; ;;     ;;  ;;     ;;    ;;          ;;;;  ;; ;;     ;;  ;;  ;;       ;;       
;;  ;;  ;; ;;;;;;;;;  ;;     ;;    ;;;;;;      ;; ;; ;; ;;     ;;  ;;   ;;;;;;  ;;;;;;   
;;  ;;  ;; ;;     ;;  ;;     ;;    ;;          ;;  ;;;; ;;     ;;  ;;        ;; ;;       
;;  ;;  ;; ;;     ;;  ;;     ;;    ;;          ;;   ;;; ;;     ;;  ;;  ;;    ;; ;;       
 ;;;  ;;;  ;;     ;; ;;;;    ;;    ;;;;;;;;    ;;    ;;  ;;;;;;;  ;;;;  ;;;;;;  ;;;;;;;; 


%% White noise analysis 
% First, I construct a nonlinearity that looks like this (it's a simple Hill function with n = 8):

stim_mean = .55;
hill_param = [1 stim_mean 8];

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
x = linspace(0,1,100);
H = hill(hill_param,x);
plot(x,H,'k')
xlabel('Stimulus')
ylabel('input nonlinearity')
prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Then, I consider a stimulus that varies from [0.2 0.9] and another that varies from [0.4 0.7], to mimic what we used in the real experiment. This is white noise, with delta correlations. I feed this stimulus through the nonlinearity, and then through a simple filter (this is my NL model). 

S = [(.4+.3*rand(1e4,1)); (.2+.7*rand(1e4,1))];
x = hill(hill_param,S);
p.A = .3; p.tau1 = 20; p.tau2 = 70; p.n = 2;
K = filter_gamma2(1:1e3,p);
R = filter(K,1,x);

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(3,1,1); hold on
plot(S)
ylabel('Stimulus')

subplot(3,1,2); hold on
plot(x)
ylabel('Output of Hill function')

subplot(3,1,3); hold on
plot(R)
xlabel('Time (s)')
ylabel('Response')
set(gca,'YLim',[0 1])

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% Now, I back out a filter from the data and use it to make linear projections, and compare the data to the linear projections. 

Khat = fitFilter2Data(S,R,'reg',1);
Rhat = filter(Khat,1,S);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(K/norm(K),'k')
plot(Khat/norm(Khat),'r')
legend({'Actual filter','Reconstructed filter'})

subplot(1,2,2); hold on
plotPieceWiseLinear(Rhat(1e3:1e4),R(1e3:1e4),'Color','b','nbins',25);
plotPieceWiseLinear(Rhat(1.1e4:end),R(1.1e4:end),'Color','r','nbins',25);
xlabel('Projected Stimulus')
ylabel('Response')


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% So it looks like it is possible for a NL model to show what we measure as variance gain control using a LN model. Note that even though there is an input nonlinearity, the output curves from the LN analysis look quite linear. 

%% White noise: effect of steepness
% How does this vary with the steepness of the input nonlinearity? In the following figure, I vary the steepness of the input nonlinearity, keeping everything else the same, and fit LN models to the data to estimate effective variance gain control as a function of input steepness. 

all_n = [1 2 4 8 16 32];

figure('outerposition',[0 0 1200 802],'PaperUnits','points','PaperSize',[1200 802]); hold on
for i = 1:length(all_n)
	% construct the NL model
	hill_param(3) = all_n(i);
	x = hill(hill_param,S);
	R = filter(K,1,x);

	% fit a LN model
	Khat = fitFilter2Data(S,R,'reg',1);
	Rhat = filter(Khat,1,S);

	subplot(2,3,i); hold on
	plotPieceWiseLinear(Rhat(1e3:1e4),R(1e3:1e4),'Color','b','nbins',25);
	plotPieceWiseLinear(Rhat(1.1e4:end),R(1.1e4:end),'Color','r','nbins',25);
	xlabel('Projected Stimulus')
	ylabel('Response')
	title(['n = ' oval(all_n(i))])

end

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% So it looks like the input nonlinearity has be steep enough for this to work. At low n, Hill functions don't seem to be able to reproduce what we see (curves with different slopes that intersect one another). 

%% White noise: effect of $k_D$
% Can we ever see something like variance gain control with $n = 1$? To check, I force $n = 1$, and vary $k_D$. 


all_k_D = logspace(log10(.1*stim_mean),log10(10*stim_mean),6);
hill_param = [1 NaN 1];

figure('outerposition',[0 0 1200 802],'PaperUnits','points','PaperSize',[1200 802]); hold on
for i = 1:length(all_k_D)
	% construct the NL model
	hill_param(2) = all_k_D(i);
	x = hill(hill_param,S);
	R = filter(K,1,x);

	% fit a LN model
	Khat = fitFilter2Data(S,R,'reg',1);
	Rhat = filter(Khat,1,S);

	subplot(2,3,i); hold on
	plotPieceWiseLinear(Rhat(1e3:1e4),R(1e3:1e4),'Color','b','nbins',25);
	plotPieceWiseLinear(Rhat(1.1e4:end),R(1.1e4:end),'Color','r','nbins',25);
	xlabel('Projected Stimulus')
	ylabel('Response')
	title(['k_D = ' oval(all_k_D(i)/stim_mean) '\mu'])

end

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% So for white noise inputs, it looks like no matter what $k_D$ is, we can never see variance gain control when $n$ is small (<8). 

 ;;;;;;   ;;;;;;;  ;;;;;;;;  ;;;;;;;;     ;;    ;;  ;;;;;;;  ;;;;  ;;;;;;  ;;;;;;;; 
;;    ;; ;;     ;; ;;     ;; ;;     ;;    ;;;   ;; ;;     ;;  ;;  ;;    ;; ;;       
;;       ;;     ;; ;;     ;; ;;     ;;    ;;;;  ;; ;;     ;;  ;;  ;;       ;;       
;;       ;;     ;; ;;;;;;;;  ;;;;;;;;     ;; ;; ;; ;;     ;;  ;;   ;;;;;;  ;;;;;;   
;;       ;;     ;; ;;   ;;   ;;   ;;      ;;  ;;;; ;;     ;;  ;;        ;; ;;       
;;    ;; ;;     ;; ;;    ;;  ;;    ;;     ;;   ;;; ;;     ;;  ;;  ;;    ;; ;;       
 ;;;;;;   ;;;;;;;  ;;     ;; ;;     ;;    ;;    ;;  ;;;;;;;  ;;;;  ;;;;;;  ;;;;;;;; 


%% Correlated Noise
% How do stimulus correlations affect this? In this section, I introuduce correlations in the white noise by a simple low-pass filter. 

stim_mean = .55;
hill_param = [1 stim_mean 8];
t_corr = 10;

S = [(.4+.3*rand(1.1e4,1)); (.2+.7*rand(1.1e4,1))];
S = filtfilt(ones(t_corr,1),t_corr,S);
S = S(1e3+1:end-1e3);
S = S - mean(S);
S = S*4;
S = S + hill_param(2);

x = hill(hill_param,S);
R = filter(K,1,x);

figure('outerposition',[0 0 1000 800],'PaperUnits','points','PaperSize',[1000 800]); hold on
subplot(3,1,1); hold on
plot(S)
ylabel('Stimulus')

subplot(3,1,2); hold on
plot(x)
ylabel('Output of Hill function')

subplot(3,1,3); hold on
plot(R)
xlabel('Time (s)')
ylabel('Response')
set(gca,'YLim',[0 1])

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now, I back out a filter from the data and use it to make linear projections, and compare the data to the linear projections. 

Khat = fitFilter2Data(S,R,'reg',1);
Rhat = filter(Khat,1,S);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(K/norm(K),'k')
plot(Khat/norm(Khat),'r')
legend({'Actual filter','Reconstructed filter'})

subplot(1,2,2); hold on
plotPieceWiseLinear(Rhat(1e3:1e4),R(1e3:1e4),'Color','b','nbins',25);
plotPieceWiseLinear(Rhat(1.1e4:end),R(1.1e4:end),'Color','r','nbins',25);
xlabel('Projected Stimulus')
ylabel('Response')


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Fits to real data
% Can we fit a NL model to real data and account for apparent variance gain control? 


[PID, LFP, fA, paradigm, orn, fly] = consolidateData(getPath(dataManager,'e30707e8e8ef6c0d832eee31eaa585aa'),1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 


% choose example
i = 8;

return

clear data
data.response = fA(global_start:global_end,i);
data.stimulus = PID(global_start:global_end,i) - min(PID(:,i));
data.response(1:1e3) = NaN;
p = fitModel2Data(@pNL,data,'nsteps',30,'p0',p);

p.Hill_n = 4;
p.Hill_K = 0.4176;
p.  tau1 = 20.3359;
p.  tau2 = 49.3125;
p.     n = 3.6250;
p.     A = 0.3324;
p.     C = 99.9383;


% clear p
% p.Hill_n = 8;
% p.Hill_K = 0.4277;
% p.  tau1 = 23.0859;
% p.  tau2 = 31.3125;
% p.     n = 3.2656;
% p.     A = 0.3334;
% p.     C = 99.9383;

% synthesize reponses using this best-fit model
S = PID(:,i) - min(PID(:,i));
R = pNL(S,p);
time = 1e-3*(1:length(S));

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,3,1:2); hold on
r2 = rsquare(R(global_start:global_end),fA(global_start:global_end,i));
plot(time,fA(:,i),'k')
plot(time,R,'r');
set(gca,'XLim',[50 180])
legend({'ab3A firing rate',['Best fit NL model, r^2 = ' oval(r2)]})
xlabel('Time (s)')
ylabel('Response (Hz)')


% back out filters
Khat = fitFilter2Data(S(global_start:global_end),R(global_start:global_end),'reg',1);
Rhat = filter(Khat,1,S);

Y = R(40e3+1:190e3);
X = Rhat(40e3+1:190e3);

Y = reshape(Y,10e3,15);
X = reshape(X,10e3,15);

subplot(1,3,3); hold on
plotPieceWiseLinear(vectorise(X(1e3:5e3,:)),vectorise(Y(1e3:5e3,:)),'Color','r','nbins',25);
plotPieceWiseLinear(vectorise(X(6e3:end,:)),vectorise(Y(6e3:end,:)),'Color','b','nbins',25);
xlabel('Projected Stimulus')
ylabel('NL model Response')


prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%% Fits to naturalistic data
% This imposes some pretty severe constraints on any NL-like model we would like to fit to the rest of our data. Can we fit an adapting NL model to the naturalistic stimulus data, with the constraint that only steep input nonlinearities are permissible?

%%
% In the following figure, I first fit a adaptive NL model where the $K_D$ can vary with the stimulus to real ab3A firing data. (a) shows the firing rate of the ab3A ORN, together with the best fit adaptive NL model. This adaptive model significantly changes its input nonlinearity (b). Then, I use this model to generate responses using the naturalistic stimulus (c). I then fit a static NL model to this synthetic data (c, blue). It fits the data well, but the static nonlinearity that I fit is much shallower (c). 

% get the naturalistic stimuli 
clear ab3 ab2
load(getPath(dataManager,'5c7dacc5b42ff0eebb980d80fec120c3'),'data','spikes')
PID = data(2).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(2).A;

% A spikes --> firing rate
fA = spiketimes2f(all_spikes,time);

tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = interp1(time,PID(i,:),tA);
end
PID = PID2; clear PID2
% some minor cleaning up
PID(end,:) = PID(end-1,:); 

PID = mean(PID,2);

NSdata.fA = mean(fA,2);
NSdata.PID = PID - min(PID);
NSdata.time = 1e-3*(1:length(NSdata.fA));

% fit an adapting NL model to this
clear data
data.response = NSdata.fA;
data.stimulus = NSdata.PID;

clear p
p.  k0 = 0.1210;
p. tau = 100;
p.   B = 1.1874;
p.tau1 = 42.1992;
p.tau2 = 400;
p.   n = 8;
p.   A = 0.2125;
p.   C = 170.6618;

% generate responses using this
[R,~,K_D] = aNLN(NSdata.PID,p);

% now fit a NL model to this
clear data
data.response = R;
data.stimulus = [NSdata.PID R];

clear p
p.k_D = 0.3702;
p.  n = 1.3086;

R2 = NLNmodel(data.stimulus,p);

time = 1e-3*(1:length(R));

figure('outerposition',[0 0 1300 600],'PaperUnits','points','PaperSize',[1300 600]); hold on
subplot(2,6,1:4); hold on
plot(time,NSdata.fA,'k')
plot(time,R,'r');
legend({'ab3A data',['adaptive NL model, r^2 = ' oval(rsquare(R,NSdata.fA))]})
set(gca,'XLim',[0 70],'YLim',[0 140])

subplot(2,6,7:10); hold on
plot(time,R,'r');
plot(time,R2,'b');
legend({'adaptive NL model',['fixed NL model, r^2 = ' oval(rsquare(R,R2))]})
set(gca,'XLim',[0 70],'YLim',[0 140])

subplot(2,6,[5 6 11 12]); hold on
clear ax
all_k_d = logspace(log10(min(K_D)),log10(max(K_D)),5);
x = logspace(-3,log10(max(NSdata.PID)),100);
for i = 1:length(all_k_d)
	ax(1) = plot(x,hill([1 all_k_d(i) 8],x),'r');
end
set(gca,'XScale','log')

ax(2) = plot(x,hill([1 p.k_D p.n],x),'b');
title('Input nonlinearities')
legend(ax,{'Actual adapting NLs','Best-fit static NL'},'Location','northwest')
xlabel('Stimulus (V)')

movePlot(gca,'right',.05)

prettyFig();

labelFigure

if being_published
	snapnow
	delete(gcf)
end


%% Version Info
%
pFooter;


