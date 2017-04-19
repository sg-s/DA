
pHeader;


% construct a DA model, vary parameters, and 

S = zeros(1e4,10);
stim_levels = logspace(-2,1,10);
for i = 1:length(stim_levels)
	S(1e3:2e3,i) = stim_levels(i);
end

clear p
p.   s0 = 0;
p.  n_z = 2;
p.tau_z = 100;
p.  n_y = 2;
p.tau_y = 10;
p.    C = 0;
p.    A = 1e3;
p.    B = 1;

S =  10 + 15*randn(1e4,1);
S = filter(ones(20,1),20,S);
S(S<0) = 0;

all_B = logspace(1.5,2,10);

noise = 0;

R = NaN*S;
for i = 1:length(all_B)
	p.B = all_B(i);
	R(:,i) = DAModelv2(S,p) + noise*randn(length(S),1);
end

% extract filters from them all
K = NaN(600,10);
for i = 1:10
	temp = fitFilter2Data(S(1e3:end),R(1e3:end,i),'filter_length',700,'offset',100);
	temp = temp(51:end-50);
	K(:,i) = temp/norm(temp);
end

% make projections and measure gain
filtertime = (1:length(K))-50;
Rp = NaN*R;
gain = NaN*all_B;
for i = 1:10
	Rp(:,i) = convolve(1:length(S),S,K(:,i),filtertime);
	ff = fit(Rp(1e3:end-1e2,i),R(1e3:end-1e2,i),'poly1');
	gain(i) = ff.p1;
end




figure('outerposition',[0 0 1100 910],'PaperUnits','points','PaperSize',[1100 910]); hold on
clear ax

% show the responses of the DA model
c = (parula(10));
ax(1) = subplot(3,4,1); hold on
for i = 1:10
	plot(R(:,i),'Color',c(i,:))
end
set(gca,'XLim',[2.5e3 3e3],'YLim',[0 40])
xlabel('Time (a.u.)')
ylabel('r(t)')

% show the backed out filters
ax(2) = subplot(3,4,2); hold on
for i = 1:10
	plot(filtertime,K(:,i),'Color',c(i,:))
end
xlabel('Lag (a.u.)')
ylabel('Filter K (norm)')

ax(3) = subplot(3,4,3); hold on
for i = 1:10
	plotPieceWiseLinear(Rp(1e3:end,i),R(1e3:end,i),'Color',c(i,:),'nbins',30);
end
xlabel('K \otimes s(t)')
ylabel('r(t)')

% show gain vs. 1/beta

ax(4) = subplot(3,4,4); hold on
for i = 1:10
	plot(1/all_B(i),gain(i),'+','Color',c(i,:))
end
l = plot(NaN,NaN,'k+');
legend(l,['r^2 = ' oval(rsquare(1./all_B,gain))],'Location','southeast');

xlabel('$1/\beta$','interpreter','latex')
ylabel('Gain (a.u.)')

;;    ;; ;;       ;;    ;;    ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;       
;;;   ;; ;;       ;;;   ;;    ;;;   ;;; ;;     ;; ;;     ;; ;;       ;;       
;;;;  ;; ;;       ;;;;  ;;    ;;;; ;;;; ;;     ;; ;;     ;; ;;       ;;       
;; ;; ;; ;;       ;; ;; ;;    ;; ;;; ;; ;;     ;; ;;     ;; ;;;;;;   ;;       
;;  ;;;; ;;       ;;  ;;;;    ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;   ;;; ;;       ;;   ;;;    ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;    ;; ;;;;;;;; ;;    ;;    ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;;;;;;; 


% OK, now we switch to a NLN model 

all_k_d = logspace(3,4,10);

% plot Hill functions
ax(5) = subplot(3,3,4); hold on
x = logspace(1,8,1e3);
for i = 1:10
	y = x./(x + all_k_d(i));
	plot(x,y,'Color',c(i,:))
end
set(gca,'XScale','log','YLim',[0 1],'XLim',[10 1e6],'XTick',[10 100 1e3 1e4 1e5]);
xlabel('Stimulus (a.u.)')
ylabel('Activity')




% generate the stimulus
S = zeros(2e4,10);
frozen_noise = randn(length(S),1);
for i = 1:10
	S(:,i) = exp(log(all_k_d(i))+frozen_noise*.1*log(all_k_d(i)));
end

% create the filter
clear p
p.n = 2;
p.tau1 = 20;
p.tau2 = 50;
p.A = .3;

K = filter_gamma2(1:500,p);
K = K/norm(K);

ax(6) = subplot(3,3,5); hold on
plot(K,'k')
xlabel('Lag (a.u.)')
ylabel('Filter K (norm)')
set(gca,'YLim',[-.02 .2])

% pass through model
R = 0*S;
a = 0*S;
for i = 1:10
	a(:,i) = S(:,i)./(S(:,i) + all_k_d(i));
	R(:,i) = filter(K,1,a(:,1));
end

ax(7) = subplot(3,3,6); hold on
for i = 1:10
	plot(R(:,i),'Color',c(i,:))
end
xlabel('Time (a.u.)')
ylabel('r(t)')
set(gca,'XLim',[2e3 2.5e3],'YLim',[3.5 4.5])

% back out filters
K = NaN(600,10);
for i = 1:10
	temp = fitFilter2Data(S(1e3:end,i),R(1e3:end,i),'filter_length',700,'offset',100,'reg',1);
	temp = temp(51:end-50);
	K(:,i) = temp/norm(temp);
end

% show the filter estimates and the actual filter
filtertime = (1:length(K))-50;
ax(8) = subplot(3,3,7); hold on
for i = 1:10
	plot(filtertime,K(:,i),'Color',c(i,:))
end
K_org = filter_gamma2(1:500,p);
K_org = K_org/norm(K_org);
plot(K_org,'k')
xlabel('Lag (a.u.)')
ylabel('Filters (norm)')

% make projections and measure gain

Rp = NaN*R;
gain = NaN*all_B;
for i = 1:10
	Rp(:,i) = convolve(1:length(S),S(:,i),K(:,i),filtertime);
	ff = fit(Rp(1e3:end-1e2,i),R(1e3:end-1e2,i),'poly1');
	gain(i) = ff.p1;
end

ax(9) = subplot(3,3,8); hold on
for i = 1:10
	plotPieceWiseLinear(Rp(1e3:end,i),R(1e3:end,i),'Color',c(i,:),'nbins',30);
end
xlabel('K \otimes s(t)')
ylabel('r(t)')

% plot gain vs. k_D
ax(10) = subplot(3,3,9); hold on
for i = 1:10
	plot(1./all_k_d(i),gain(i),'+','Color',c(i,:))
end
l = plot(NaN,NaN,'k+');
legend(l,['r^2 = ' oval(rsquare(1./all_k_d,gain))],'Location','southeast')


xlabel('$1/K_{D}$','interpreter','latex')
ylabel('Gain (a.u.)')

prettyFig('fs',15)

ax(5).Position(2) = .38;
ax(6).Position(2) = .38;
ax(7).Position(2) = .38;

ax(8).Position(2) = .08;
ax(9).Position(2) = .08;
ax(10).Position(2) = .08;

labelFigure('x_offset',0,'y_offset',.01)

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


