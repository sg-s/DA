% Fast Gain Control
% 
% created by Srinivas Gorur-Shandilya at 3:28 , 24 September 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

pHeader;


%         ########    ###     ######  ########     ######      ###    #### ##    ## 
%         ##         ## ##   ##    ##    ##       ##    ##    ## ##    ##  ###   ## 
%         ##        ##   ##  ##          ##       ##         ##   ##   ##  ####  ## 
%         ######   ##     ##  ######     ##       ##   #### ##     ##  ##  ## ## ## 
%         ##       #########       ##    ##       ##    ##  #########  ##  ##  #### 
%         ##       ##     ## ##    ##    ##       ##    ##  ##     ##  ##  ##   ### 
%         ##       ##     ##  ######     ##        ######   ##     ## #### ##    ## 

%          ######   #######  ##    ## ######## ########   #######  ##       
%         ##    ## ##     ## ###   ##    ##    ##     ## ##     ## ##       
%         ##       ##     ## ####  ##    ##    ##     ## ##     ## ##       
%         ##       ##     ## ## ## ##    ##    ########  ##     ## ##       
%         ##       ##     ## ##  ####    ##    ##   ##   ##     ## ##       
%         ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##       
%          ######   #######  ##    ##    ##    ##     ##  #######  ######## 


%% Fast Gain Control in ORNs
% 

fig_handle=figure('outerposition',[0 0 900 1000],'PaperUnits','points','PaperSize',[900 1000]); hold on
clf(fig_handle);
axes_handles(1) = subplot(9,3,[1 2 4 5]); 
axes_handles(2) = subplot(9,3,[7 8 10 11]); 
axes_handles(3) = subplot(9,3,[13 14 16 17]);
axes_handles(4) = subplot(9,3,[19 22 25]);
axes_handles(5) = subplot(9,3,[20 23 26]);
axes_handles(6) = subplot(9,3,[21 24 27]);

for i = 1:length(axes_handles)
	hold(axes_handles(i),'on');
end

movePlot(axes_handles(1),'up',.03)
movePlot(axes_handles(2),'up',.04)
movePlot(axes_handles(3),'up',.02)

movePlot(axes_handles(4),'down',.02)
movePlot(axes_handles(5),'down',.02)
movePlot(axes_handles(6),'down',.02)

% make three more axes for explanation of how we calculate instantaneous gain
inst_ax(1) = subplot(18,3,6); hold on
inst_ax(2) = subplot(18,3,12); hold on
inst_ax(3) = subplot(9,3,12); hold on

% resize some plots and move them into position
temp = get(inst_ax(3),'Position');
temp(4) = .2;
set(inst_ax(3),'Position',temp);

for i = 1:2
	temp = get(inst_ax(i),'Position');
	temp(4) = 2*temp(4);
	set(inst_ax(i),'Position',temp);
end

movePlot(inst_ax(2),'down',.3)

% make all plots smaller
for i = 1:3
	temp = get(inst_ax(i),'Position');
	temp(3) = .18;
	set(inst_ax(i),'Position',temp);
	movePlot(inst_ax(i),'right',.05)
end


load('/local-data/DA-paper/large-variance-flicker/ab3/2015_01_28_CS_ab3_2_EA.mat')
PID = data(4).PID;
time = 1e-4*(1:length(PID));
all_spikes = spikes(4).A;
B_spikes = spikes(4).B;
load('/local-data/DA-paper/large-variance-flicker/ab3/2015_01_28_CS_ab3_3_EA.mat')
PID = vertcat(PID,data(4).PID);
all_spikes = vertcat(all_spikes,spikes(4).A);

% A spikes --> firing rate
hash = dataHash(full(all_spikes));
cached_data = cache(hash);
if isempty(cached_data)
	disp('Computing firing rate for A neuron...')
	fA = spiketimes2f(all_spikes,time);
	cache(hash,fA);
else
	fA = cached_data;
end


tA = 1e-3*(1:length(fA));
PID2 = fA;
for i = 1:width(PID2)
	PID2(:,i) = PID(i,1:10:end);
end
PID = PID2; clear PID2


% set up a colour map
c = parula(8);

% plot stimulus
plot(axes_handles(1),tA,mean(PID,2),'k')
set(axes_handles(1),'XLim',[0 60],'XTick',[])
ylabel(axes_handles(1),'Stimulus (V)')

% make linear predictions for the whole experiment
[this_K, filtertime_full] = fitFilter2Data(mean(PID(10e3:end,:),2),mean(fA(10e3:end,:),2),'reg',1,'filter_length',1200,'offset',200);
K = this_K(100:end-101);
filtertime = filtertime_full(100:end-101);

fp = NaN*fA;
for i = 1:width(fA)
	fp(:,i) = convolve(tA,PID(:,i),K,filtertime);
end

clear l
R = mean(fA,2);
[ax,plot1,plot2] = plotyy(axes_handles(2),tA,R,tA,mean(fp,2));
set(ax(1),'XLim',[0 60],'YLim',[min(mean(fA,2)) 1.3*max(mean(fA,2))])
set(ax(2),'XLim',[0 60],'YLim',[min(mean(fp,2)) max(mean(fp,2))])
set(ax(2),'YTick',-.1:.1:.5,'YColor','r')
set(plot1,'Color','k')
set(plot2,'Color','r')
set(ax(1),'YColor',[0 0 0])
ylabel(ax(1),'ORN Response (Hz)')
ylabel(ax(2),'Projected Stimulus (V)')
set(axes_handles(2),'box','off')



hl_min = 1e-2;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);

resp = mean(fA(10e3:55e3,[3:10 13:20]),2);
pred = mean(fp(10e3:55e3,[3:10 13:20]),2);
time = 1e-3*(1:length(resp));
stim = mean(PID(10e3:55e3,[3:10 13:20]),2);


% plot gain vs preceding stimulus
[mean_stim,inst_gain,e] = findInstGain(mean(PID,2),mean(fA,2),mean(fp,2),50);
rm_this = (isnan(mean_stim) | isnan(inst_gain)) | e < .8 | inst_gain < 0;
mean_stim(rm_this) = NaN;
inst_gain(rm_this) = NaN;
mean_stim(1:1e3) = NaN; inst_gain(1:1e3) = NaN;

% explain how we compute the inst. gain
plot(inst_ax(1),tA,mean(fA,2),'k')
plot(inst_ax(2),tA,mean(fp,2),'r')
set(inst_ax(1),'XLim',[21 22],'XTick',[],'XTickLabel',{},'XColor','w','YLim',[0 60])
set(inst_ax(2),'XLim',[21 22],'XTick',[],'XTickLabel',{},'XColor','w','YLim',[0 .4])
% insert a marker at 21.5
plot(inst_ax(1),[21.5 21.5],[0 200],'k')
plot(inst_ax(2),[21.5 21.5],[0 0.5],'k')

x = mean(fp,2); y = mean(fA,2);
x = x(21.5e3-50:21.5e3); y = y(21.5e3-50:21.5e3);
plot(inst_ax(3),x,y,'b.')
ff = fit(x(:),y(:),'poly1');
plot(inst_ax(3),sort(x),ff(sort(x)),'b')

plot(inst_ax(1),tA(21.5e3-50:21.5e3),y,'k','LineWidth',4)
plot(inst_ax(2),tA(21.5e3-50:21.5e3),x,'r','LineWidth',4)

set(inst_ax(3),'XColor','r')
set(inst_ax(2),'YColor','r')

temp = get(inst_ax(3),'Position');
set(inst_ax(3),'Position',[.74 .59 .17 .18])



plot(axes_handles(3),tA(1:10:end),inst_gain(1:10:end),'b.')
ylabel(axes_handles(3),'Inst. Gain (Hz/V)')
xlabel(axes_handles(3),'Time (s)')
set(axes_handles(3),'XLim',[0 60],'YScale','log','YLim',[10 1e4])

% show that the inst. gain varies a bit
[y,x] = hist(log(inst_gain),100);
y = y/sum(y);
plot(axes_handles(4),exp(x),y,'b')
set(axes_handles(4),'XScale','log','XLim',[10 1e4],'XTick',[1e1 1e2 1e3 1e4],'YLim',[0 max(y)*1.1])
xlabel(axes_handles(4),'Instantaneous Gain (Hz/V)')
ylabel(axes_handles(4),'Probability')


% show the inst gain vs. the mean stim
K = ones(500,1);
% filter the stimulus
mean_stim = abs(filter(K,sum(K),mean(PID,2)));
mean_stim(1:1e3) = NaN;
[~,data] = plotPieceWiseLinear(mean_stim,inst_gain,'nbins',100,'Color','b','use_std',false,'make_plot',false);
l = errorShade(axes_handles(5),data.x,data.y,data.ye,'Color',[0 0 1]);
temp1 = mean_stim(~isnan(inst_gain));
temp2 = inst_gain(~isnan(inst_gain));
legend(l(1),['\rho =' oval(spear(temp1(1:10:end),temp2(1:10:end)),3)])
xlabel(axes_handles(5),'Stimulus in preceding 500ms (V)')
ylabel(axes_handles(5),'Inst. Gain (Hz/V)')
set(axes_handles(5),'YScale','log','XLim',[min(mean_stim) max(mean_stim)],'YLim',[1e1 1e4])



% now vary the gain filter and show it is the best possible
rho = NaN*history_lengths;
for i = 1:length(history_lengths)
	hl = floor(history_lengths(i)*1e3);
	K = ones(hl,1);
	% filter the stimulus
	shat = abs(filter(K,sum(K),mean(PID,2)));

	% remove some junk
	temp = inst_gain(~isnan(inst_gain));
	shat = shat(~isnan(inst_gain));

	% use the Spearman rank correlation
	rho(i) = (spear(shat(1:10:end),temp(1:10:end)));
end

plot(axes_handles(6),history_lengths,rho,'b+')
set(axes_handles(6),'XScale','log','XMinorTick','on','XLim',[1e-2 1e1])
ylabel(axes_handles(6),'\rho')
xlabel(axes_handles(6),'History Length (s)')

prettyFig('plw=1.5;','lw=1.2;','fs=14;','FixLogX=true;')

ylabel(inst_ax(1),'ORN Response (Hz)','FontSize',10)
ylabel(inst_ax(2),'Proj. Stimulus (V)','FontSize',10)
set(inst_ax(1),'FontSize',10)
set(inst_ax(2),'FontSize',10)

if being_published
	snapnow
	delete(gcf)
end


%% Supplementary Figure
% 

figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on

[mean_stim,inst_gain,e] = findInstGain(mean(PID,2),mean(fA,2),mean(fp,2),50);
subplot(2,3,1), hold on
[y,x] = hist(e,50);
y = y/sum(y);
plot(x,y,'r')
plot([.8 .8],[0 .5],'k--')
xlabel('r^2')
ylabel('p(r^2)')

subplot(2,3,2), hold on
plot([min(mean_stim) max(mean_stim)],[.8 .8],'k--')
l = plot(mean_stim(1:50:end),e(1:50:end),'r.');
xlabel('Stimulus in preceding 50ms (V)')
legend(l,['r^2 = ' oval(rsquare(mean_stim,e))],'Location','southeast')
ylabel('r^2')

subplot(2,3,3), hold on
plot([min(inst_gain) max(inst_gain)],[.8 .8],'k--')
l = plot(inst_gain(1:50:end),e(1:50:end),'r.');
set(gca,'XLim',[-5e3 1e4])
xlabel('Inst. gain (Hz/V)')
legend(l,['r^2 = ' oval(rsquare(inst_gain,e))],'Location','southeast')
ylabel('r^2')


clear ph
ph(3) = subplot(2,3,4); hold on
ph(4) = subplot(2,3,5); hold on


hl_min = 1e-1;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);

[~,~,~,~,~,history_lengths]=gainAnalysisWrapper('response',resp,'prediction',pred,'stimulus',stim,'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',.5,'use_cache',true,'engine',@gainAnalysis);
set(ph(4),'XLim',[.09 11]) % to show .1 and 10 on the log scale
set(ph(3),'XLim',[min(pred) max(pred)])
xlabel(ph(3),'Projected Stimulus (V)')
ylabel(ph(3),'ORN Response (Hz)')

% show the p-value
axes(ph(3))
text(.1,60,'p < 0.01')
title(ph(3),[])

% fix some labels
set(ph(4),'YScale','log','YTick',[0.5 1 2],'YLim',[0.4 2.5],'XLim',[0.09 10.1],'YMinorTick','on')
title(ph(4),'Linear Model')


% now do the same analysis with the LN model

fp = mean(fp,2);
temp = fit(fp(~(isnan(fp) | isnan(R))),R(~(isnan(fp) | isnan(R))),'poly1');
fp = fp*temp.p1;
fp = fp+temp.p2;

clear p
p.A = 57.1534;
p.k = 23.6690;
p.n = 2.9341;
fp_hill = hill(p,fp);

clear ph
ph(4) = subplot(2,3,6); hold on

hl_min = .1;
hl_max = 10;
history_lengths = [logspace(log10(hl_min),log10(.5),15) logspace(log10(.5),log10(10),15)];
history_lengths = unique(history_lengths);

resp = mean(fA(10e3:55e3,[3:10 13:20]),2);
pred = (fp_hill(10e3:55e3));
time = 1e-3*(1:length(resp));
stim = mean(PID(10e3:55e3,[3:10 13:20]),2);

[p,~,~,~,~,history_lengths] = gainAnalysisWrapper('response',resp,'prediction',pred,'stimulus',stim,'time',time,'ph',ph,'history_lengths',history_lengths,'example_history_length',.5,'use_cache',true,'engine',@gainAnalysis);
set(ph(4),'XLim',[.09 11]) % to show .1 and 10 on the log scale
title(ph(4),'LN Model')

prettyFig('FixLogX=1;','fs=18;')

if being_published
	snapnow
	delete(gcf)
end

%% Gain Control beyond a static nonlinearity
% In this section, we fit a static nonlinearity to the data, and then repeat the inst. gain analysis to show that there is still some gain control that cannot be explained by a static output nonlinearity. 

load('/Users/sigbhu/code/da/data/LVP_orn_data.mat','od')
figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end
hl = [50 300 1e4];
for i = 1:3
	plot_handles = plot(od,[ax(i) ax(4)],'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',hl(i),'nbins',300);
	title(ax(i),['\tau_{H} = ' oval(hl(i)) 'ms'])
	if i < 3
		delete(plot_handles(2).f2)
	end
	ylabel(ax(i),'Inst Gain (norm)')
end
set(ax(1:3),'XScale','log','YScale','log','YLim',[.5 10])
set(ax(4),'YLim',[-1 1])

prettyFig('fs=14;','FixLogY=true;');

%%
% As is clear from the figure above, fitting an output nonlinearity eliminates most of the gain control that we see. We did fit a nonlinearity when we used the Clark et al. style fast gain analysis, and we still found significant dynamic gain control. Since this analysis method is  less sensitive (because we don't bin data), it doesn't see a lot of gain control after fitting a nonlinearity.  



%% Version Info
% 

pFooter;
