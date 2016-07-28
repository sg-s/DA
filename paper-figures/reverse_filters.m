

pHeader;


clear cdata
cdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
cdata = cleanMSGdata(cdata);

v2struct(cdata)

time = 1e-3*(1:length(PID));

%% Predictive Filters
% In this document, I attempt to back out filters that predict the stimulus given the response. This obviously will not work if the spectral content of the input and the output are very different: for example, white noise inputs can't be predicted from low-passed versions of the white noise, no matter what. To consider a realistic case, I consider the correlated fluctuating stimulus I used to stimualate ab3A ORNs:

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
S = PID(:,11);
plot(time,S);

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now, I consider an idealized filter and generate some responses using this. I then flip the time series around and try to back out a filter from the response to the stimulus. 

clear p
p.tau2 = 70;
p.tau1 = 20;
p.   n = 2;
p.   A = 0.3000;
K = filter_gamma2(1:500,p);
R = filter(K,1,S);

a = 15e3;
z = 45e3;
K_rev = fitFilter2Data(flipud(R(a:z)),flipud(S(a:z)),'reg',1,'offset',200,'filter_length',1000);
K_rev(1:100) = [];
filtertime = 1e-3*(1:length(K_rev)) - .1;
sp = convolve(time,flipud(R),K_rev,filtertime);

figure('outerposition',[0 0 1501 500],'PaperUnits','points','PaperSize',[1501 500]); hold on
subplot(1,3,1); hold on
plot(K)
xlabel('Lag (ms)')
title('Forward filter')

subplot(1,3,2)
plot(filtertime,K_rev)
xlabel('Anti-lag (s)')
title('Reverse filter')

subplot(1,3,3); hold on
plot(flipud(sp),S,'k.')
title(['r^2 = ' oval(rsquare(flipud(sp),S))])
xlabel('K_{rev} \otimes R')
ylabel('S')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% That seems to work. Now, we apply it to the real data. 

K_rev = NaN(800,length(paradigm));

for i = 1:width(PID)
	if paradigm(i) == 1
		S = PID(:,i);
		R = fA(:,i);
		K_rev(:,i) = fitFilter2Data(flipud(R(a:z)),flipud(S(a:z)),'reg',1,'offset',200,'filter_length',800);
	end
end

K_rev(1:100,:) = [];
filtertime = 1e-3*(1:length(K_rev)) - .1;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(filtertime,nanmean(K_rev,2),'r')
plot(filtertime,nanmean(K2(:,paradigm==1),2),'b')
legend({'K_{rev}','K_{forward}'})
xlabel('Lag (or reverse lag) (s)')
ylabel('Filter amplitude')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% As predicted, the reverse filters are simply the forward filters flipped around. 




%% Version Info
%
pFooter;

