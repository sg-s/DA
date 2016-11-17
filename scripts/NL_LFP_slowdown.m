
pHeader;

%% Can NL models reproduce apparent LFP slowdown?
% In this document, I attempt to first fit the LFP data using a NL model, and then perform the slowdown analysis on the LFP, to see if the LFP does indeed slowdown. The following figure shows the LFP and the best-fit NL model. 

% get the data
load('/local-data/DA-paper/data-for-paper/fig7/nat-stim-ab3/2016_02_09_CS_F2_ab3_5.mat')
LFP = data(2).voltage(:,1:10:end)';
LFP = mean(LFP,2);
LFP = LFP - mean(LFP(1:5e3));
PID = data(2).PID(:,1:10:end)';
PID = mean(PID,2);
PID = PID - min(PID(1:5e3));

% fit the NL model to it
clear data p
data.response = LFP;
data.stimulus = [PID LFP];
p.k_D = 0.0835;
p.n = 1;

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
time = 1e-3*(1:length(LFP));
plot(time,data.response,'k')
[R,K] = NLModel(data.stimulus,p);
plot(time,R,'r')
legend({'Data',['NL model, r^2 =' oval(rsquare(data.response,R))]},'Location','southeast')
xlabel('Time (s)')
ylabel('\DeltaLFP (mV)')
set(gca,'XLim',[0 70])

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%%
% Now we perform the analysis of whether the LFP slows down or not, with both the real LFP and the NL model generated LFP. It looks like somehow, the NL model can reproduce what we see as a slowdown in the LFP. This means this is probably not a real effect. In addition, I plot the stimulus autocorrelation time (measured in 1 s windows), where the stimulus autocorrelation time is defined as the lag when the stimulus autocorrelation function drops below 1/e. 


min_acceptable_corr = .5;
min_acceptable_lag = 2;


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on

[lag, mean_x, max_corr] = findLagAndMeanInWindow(PID(5e3:end-5e3),-LFP(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color','k','nbins',19);

[lag, mean_x, max_corr] = findLagAndMeanInWindow(PID(5e3:end-5e3),-R(5e3:end-5e3),1e3,50);
rm_this = lag<min_acceptable_lag | max_corr < min_acceptable_corr;
lag(rm_this) = [];
mean_x(rm_this) = [];
plotPieceWiseLinear(mean_x,lag,'Color','r','nbins',19);

% find auto correlation time
act = NaN*PID;
mean_x = NaN*PID;
for i = 5e3:50:(length(PID)-5e3)
	x = PID(i-5e2:i+5e2);
	a = autocorr(x,length(x)-1);
	act(i) = find(a<1/exp(1),1,'first');
	mean_x(i) = mean(x);
end
rm_this = isnan(mean_x) | isnan(act);
mean_x(rm_this) = []; act(rm_this) = [];
plotPieceWiseLinear(mean_x,act,'Color','g','nbins',19);

xlabel('\mu_{Stimulus} in preceding 1s (V)')
ylabel('Lag (ms)')

legend({'Data','NL model','Stimulus auto corr.'},'Location','northwest')

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;



