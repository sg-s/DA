
pHeader;

% make a filter, stimulus and some response 
S = randn(1e4,1);
clear p
p.A = .7;
p.tau1 = 10;
p.tau2 = 30;
p.n = 1;
K = filter_gamma2(1:300,p);
R = filter(K,1,S + .2*randn(1e4,1));

figure('outerposition',[0 0 1200 701],'PaperUnits','points','PaperSize',[1200 701]); hold on
subplot(2,3,1); hold on
plot(S,'k')
xlabel('t')
ylabel('s(t)')
set(gca,'XLim',[1e3 2e3])
subplot(2,3,4); hold on
plot(R,'r')
xlabel('t')
ylabel('r(t)')
set(gca,'XLim',[1e3 2e3])


% plot distributions 
subplot(2,3,2); hold on
[hy,hx] = histcounts(S,100); 
hy = hy/sum(hy); hy = hy/mean(diff(hx));
hx = hx(1:end-1) + mean(diff(hx));
plot(hx,hy,'k')
[hy,hx] = histcounts(R,100);
hy = hy/sum(hy); hy = hy/mean(diff(hx));
hx = hx(1:end-1) + mean(diff(hx)); 
plot(hx,hy,'r')
xlabel('Value')
ylabel('p.d.f')
legend({'s','r'})

% plot auto correlation functions 
subplot(2,3,5); hold on
[y,lags] = autocorr(S,1e3);
plot(lags,y,'k')
[y,lags] = autocorr(R,1e3);
plot(lags,y,'r')
set(gca,'XScale','log')
xlabel('Lag')
ylabel('Autocorrelation')
legend({'s(t)','r(t)'})

% plot actual filter, recovered filter, cross correlation 
subplot(2,3,6); hold on

R = reshape(R,5e2,20);
S = reshape(S,5e2,20);
x = NaN(1e3-1,20);
for i = 1:20
	x(:,i) = xcorr(R(:,i),S(:,i));
end

x = mean(x,2);
filtertime = (0:length(x)-1) - 5e2;
plot(filtertime,x/norm(x),'g')

Khat = fitFilter2Data(S(:),R(:),'reg',0,'offset',100,'filter_length',400);
filtertime = (0:length(Khat)-1) - 100;
plot(filtertime,Khat/norm(Khat),'b')

filtertime = 0:length(K)-1;
plot(filtertime,K/norm(K),'r')

set(gca,'XLim',[-100 500])
xlabel('Lag')
ylabel('Filter amplitude (norm)')
legend({'Cross correlation','estimated filter','Actual filter'})


% show the covariance matrix 
% chop up the stimulus into blocks  
S = S(:);
filter_length = 300;
s = zeros(length(S(:))-filter_length+1, filter_length);
only_these_points = filter_length+1:length(S);
for i = 1:length(only_these_points)
	s(i,:) = S(only_these_points(i):-1:only_these_points(i)-filter_length+1);
end

% compute covariance matrix
C = s'*s; % this is the covariance matrix, scaled by the size of the C
% scale reg by mean of eigenvalues

subplot(2,3,3); hold on
imagesc(C)
axis square
colorbar
title('C')
set(gca,'XLim',[1 30],'YLim',[1 30])
axis off
axis ij

prettyFig()

labelFigure('column_first',true,'x_offset',0);

if being_published	
	snapnow	
	delete(gcf)
end


% NOW show the same thing with correlated stimulus 


S = randn(5e3,1);
clear p
p.A = .7;
p.tau1 = 10;
p.tau2 = 30;
p.n = 1;
K = filter_gamma2(1:300,p);
R = filter(K,1,S + .2*randn(length(S),1));

act = 50;
S2 = filter(ones(30,1),act,S);
R2 = filter(K,1,S2 + .2*randn(length(S),1));

figure('outerposition',[0 0 1300 900],'PaperUnits','points','PaperSize',[1300 900]); hold on
subplot(3,4,1); hold on
plot(S,'k')
xlabel('t')
ylabel('s(t)')
set(gca,'XLim',[1e3 2e3])

subplot(3,4,5); hold on
plot(S2,'k')
xlabel('t')
ylabel('s_2(t)')
set(gca,'XLim',[1e3 2e3])


% plot auto correlation functions 
subplot(3,4,2); hold on
[y,lags] = autocorr(S,1e3);
plot(lags,y,'k')
set(gca,'XScale','log')
xlabel('Lag')
ylabel('Autocorrelation')
set(gca,'YLim',[-.2 1.1])


subplot(3,4,6); hold on
[y,lags] = autocorr(S2,1e3);
plot(lags,y,'k')
set(gca,'XScale','log')
xlabel('Lag')
ylabel('Autocorrelation')
set(gca,'YLim',[-.2 1.1])

% show the covariance matrices 

filter_length = 300;
s = zeros(length(S(:))-filter_length+1, filter_length);
only_these_points = filter_length+1:length(S);
for i = 1:length(only_these_points)
	s(i,:) = S(only_these_points(i):-1:only_these_points(i)-filter_length+1);
end
s2 = zeros(length(S2(:))-filter_length+1, filter_length);
only_these_points = filter_length+1:length(S);
for i = 1:length(only_these_points)
	s2(i,:) = S2(only_these_points(i):-1:only_these_points(i)-filter_length+1);
end

C = s'*s; 
C2 = s2'*s2; 
C3 = C2 + mean(eig(C2))*eye(length(C));

subplot(3,4,3); hold on
imagesc(C)
axis square
colorbar
title('$C$','interpreter','latex')
set(gca,'XLim',[1 30],'YLim',[1 30])
axis off
axis ij

subplot(3,4,7); hold on
imagesc(C2)
axis square
colorbar
title('$C_2$','interpreter','latex')
set(gca,'XLim',[1 30],'YLim',[1 30])
axis off
axis ij

subplot(3,4,11); hold on
imagesc(C3)
axis square
colorbar
title('$C_2 + rI$','interpreter','latex')
set(gca,'XLim',[1 30],'YLim',[1 30])
axis off
axis ij

% show the filters 
K1 = fitFilter2Data(S,R,'reg',0,'filter_length',300);
K2 = fitFilter2Data(S2,R2,'reg',0,'filter_length',300);
K3 = fitFilter2Data(S2,R2,'reg',1,'filter_length',300);

subplot(3,4,4); hold on
filtertime = 0:length(K)-1;
plot(filtertime,K/norm(K),'k')
plot(filtertime,K1/norm(K1),'r')
xlabel('Filter lag')
ylabel('Filter (norm)')
legend({'Actual','estimated'})

subplot(3,4,8); hold on
filtertime = 0:length(K)-1;
plot(filtertime,K/norm(K),'k')
plot(filtertime,K2/norm(K2),'r')
xlabel('Filter lag')
ylabel('Filter (norm)')
legend({'Actual','estimated'})

subplot(3,4,12); hold on
filtertime = 0:length(K)-1;
plot(filtertime,K/norm(K),'k')
plot(filtertime,K3/norm(K3),'r')
xlabel('Filter lag')
ylabel('Filter (norm)')
legend({'Actual','estimated'})

prettyFig();

labelFigure('column_first',false,'x_offset',0);

subplot(3,4,12);
box off

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;


