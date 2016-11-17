
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
% Now we perform the analysis of whether the LFP slows down or not, with both the real LFP and the NL model generated LFP. It looks like somehow, the NL model can reproduce what we see as a slowdown in the LFP. This means this is probably not a real effect. 


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

xlabel('\mu_{Stimulus} in preceding 1s (V)')
ylabel('Lag (ms)')

legend({'Data','NL model'})

prettyFig()

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;



