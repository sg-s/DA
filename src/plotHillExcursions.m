% plotHillExcursions.m
% 
% created by Srinivas Gorur-Shandilya at 2:12 , 04 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function plotHillExcursions(S,x,y,fix_n)


% normalise 
x = x - min(x); x = x/max(x);
y = y - min(y); y = y/max(y);

[ons,offs] = computeOnsOffs(y>.1);

r2 = 0*ons; k = 0*ons; n = 0*ons;

% first fit to all the data
ft = fittype('hill2(x,k,n)');
ff = fit(x(~isnan(x)),y(~isnan(x)),ft,'StartPoint',[.5 1],'Upper',[1 100],'Lower',[0 1]);
n0 = ff.n;
k0 = ff.k;
for i = 1:length(ons)
	temp1 = x(ons(i):offs(i));
	temp2 = y(ons(i):offs(i));

	if fix_n
		ff = fit(temp1,temp2,ft,'StartPoint',[k0 n0],'Lower',[0 n0],'Upper',[1 n0],'MaxIter',1e3);
	else
		ff = fit(temp1,temp2,ft,'StartPoint',[k0 n0],'Lower',[0 1],'Upper',[1 100],'MaxIter',1e3);
	end
	r2(i) = rsquare(temp2,ff(temp1));
	n(i) = ff.n;
	k(i) = ff.k;
end

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:5
	ax(i) = subplot(2,3,i); hold on
end
c = lines(length(ons));
for i = 1:length(ons)
	if r2(i) > .5
		temp1 = x(ons(i):offs(i));
		temp2 = y(ons(i):offs(i));
		plot(ax(1),temp1,temp2,'.','Color',c(i,:));
		plot(ax(2),[0:.01:1],hill2([0:.01:1],k(i),n(i)),'Color',c(i,:));
	end
end
set(ax(1),'YLim',[0 1])
ylabel(ax(1),'ORN Response (norm)')
xlabel(ax(1),'Proj. Stim (norm)')

xlabel(ax(2),'Proj. Stimulus (norm)')
ylabel(ax(2),'Hill function output (norm)')


plot(ax(3),r2,k,'k+')
set(ax(3),'XLim',[0 1],'YLim',[0 .5])
xlabel(ax(3),'r^2')
ylabel(ax(3),'K_D')

% find stimulus in preceding 500ms
shat = computeSmoothedStimulus(S,500);
sx = 0*k;
for i = 1:length(ons)
	sx(i) = mean(shat(ons(i):offs(i)));
end
l = plot(ax(4),sx(r2>.5),k(r2>.5),'k+');
legend(l,['\rho = ' oval(spear(sx(r2>.5),k(r2>.5)))],'Location','southeast')
xlabel(ax(4),'\mu_{stimulus} in preceding 500ms')
ylabel(ax(4),'K_D')
set(ax(4),'YLim',[0 .5])

history_lengths = round(logspace(1,4,50));
for j = 1:length(history_lengths)
	shat = computeSmoothedStimulus(S,history_lengths(j));
	sx = 0*k;
	for i = 1:length(ons)
		sx(i) = mean(shat(ons(i):offs(i)));
	end
	s = spear(sx(r2>.5),k(r2>.5));
	plot(ax(5),history_lengths(j),s,'k+')
end
set(ax(5),'XScale','log')
xlabel(ax(5),'History Length (ms)')
ylabel(ax(5),'\rho')
