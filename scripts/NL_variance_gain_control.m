
pHeader;

%% Can NL models typically lead to variance gain control?
% In this document, I check if a simple NL model can lead to something that we think of as variance gain control (along the lines of Nemenman's paper). 

%%
% First, I construct a nonlinearity that looks like this (it's a simple Hill function with n = 8):

hill_param = [1 .55 8];

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
% Then, I consider a stimulus that varies from [0.2 0.9] and another that varies from [0.4 0.7], to mimic what we used in the real experiment. I feed this stimulus through the nonlinearity, and then through a simple filter. 

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

%%
% So it looks like it is possible for a NL model to show what we measure as variance gain control using a LN model. How does this vary with the steepness of the input nonlinearity? In the following figure, I vary the steepness of the input nonlinearity, keeping everything else the same, and fit LN models to the data to estimate effective variance gain control as a function of input steepness. 

all_n = [1 2 4 8 16 32];

figure('outerposition',[0 0 1200 802],'PaperUnits','points','PaperSize',[1200 802]); hold on
for i = 1:length(all_n)
	% construct the NL model
	hill_param(3) = all_n(i);
	x = hill(hill_param,S);
	R = filter(K,1,x);

	% fit a LN model
	Khat = fitFilter2Data(S,R,'reg',0);
	Rhat = filter(Khat,1,S);

	subplot(2,3,i); hold on
	plotPieceWiseLinear(Rhat(1e3:5e3),R(1e3:5e3),'Color','b','nbins',25);
	plotPieceWiseLinear(Rhat(6e3:end),R(6e3:end),'Color','r','nbins',25);
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



%% Version Info
%
pFooter;


