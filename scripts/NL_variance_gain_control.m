
pHeader;

%% Can NL models typically lead to variance gain control?
% In this document, I check if a simple NL model can lead to something that we think of as variance gain control (along the lines of Nemenman's paper). 

%%
% First, I construct a nonlinearity that looks like this (it's a simple Hill function):

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
x = linspace(0,10,100);
H = hill([1 5 8],x);
plot(x,H,'k')
xlabel('Stimulus')
ylabel('input nonlinearity')
prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Then, I consider a stimulus that varies from [4 6] and another that varies from [2 8]. I feed this stimulus through the nonlinearity, and then through a simple filter. 

S = [(4+2*rand(1e4,1)); (2+6*rand(1e4,1))];
x = hill([1 5 8],S);
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
set(gca,'YLim',[.2 .5])

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end

%%
% Now, I back out a filter from the data and use it to make linear projections, and compare the data to the linear projections. 

Khat = fitFilter2Data(S,R,'reg',0);
Rhat = filter(Khat,1,S);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(K/norm(K),'k')
plot(Khat/norm(Khat),'r')
legend({'Actual filter','Reconstructed filter'})

subplot(1,2,2); hold on
plotPieceWiseLinear(Rhat(1e3:5e3),R(1e3:5e3),'Color','b','nbins',25);
plotPieceWiseLinear(Rhat(6e3:end),R(6e3:end),'Color','r','nbins',25);
xlabel('Projected Stimulus')
ylabel('Response')


prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


