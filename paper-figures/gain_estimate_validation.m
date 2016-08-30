pHeader

%% How to estimate gain 
% In this document, we use the DA model to validate various ways of estimating gain, and confirming the absence or presence of Weber's Law. 

clear cdata
cdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
cdata = cleanMSGdata(cdata);


% load DA model 
% generate neuron-specific DA model predictions
load('DA_model_fit_to_MSG.mat','p')
cdata2 = cdata;
for i = 1:width(cdata.PID)
	p(i).s0 = 0;
	cdata2.fA(:,i) = DAModelv2(cdata.PID(:,i),p(cdata.orn(i)));
end

cdata2.LFP = NaN*cdata.LFP;

% extract filters, etc
cdata2 = cleanMSGdata(cdata2);


mean_stim = mean(cdata.PID(35e3:55e3,:));


% also generate responses to simple steps 
cdata3 = cdata2;
back_stim = logspace(log10(min(mean_stim)),log10(max(mean_stim)),size(cdata.PID,2));
for i = 1:length(cdata.paradigm)
	S = ones(60e3,1)*back_stim(i);
	S(30e3:31e3) = 2*S(30e3:31e3);
	cdata3.PID(:,i) = S;
	R = DAModelv2(S,p(1));

	% compute gain
	dR = (max(R(30e3:31e3))-mean(R(28e3:29e3)));
	dS = back_stim(i);
	cdata3.fA_gain(i) = dR/dS;

end

%% Initial tests
% In this section, we use the DA model operating on experimental white noise stimuli with increasing means to generate synthetic responses. We then push this synthetic data through the normal analysis pipeline, where we extract filters, estimate gain, etc. The result of this is shown on the left panel. In comparison, we compute gains from data where the DA model responses to pulses on backgrounds. Here, gain is computed from the fractional response, and no filter estimation is needed. In both cases, we compare this to the Weber prediction (red line).  

% plot gain of DA model vs mean stim
figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
plot(mean_stim,cdata2.fA_gain,'k+')
set(gca,'XScale','log','YScale','log','XLim',[.15 1.5],'YLim',[2 200])
cf = fit(vectorise(mean_stim),cdata2.fA_gain,'power1','Upper',[Inf -1],'Lower',[-Inf -1]);
plot(mean_stim,cf(mean_stim),'r')
xlabel('Mean (white noise stimulus)')
ylabel('Gain (Hz/V)')
title('Gain estimation using white noise')

subplot(1,2,2); hold on
plot(back_stim,cdata3.fA_gain,'k+')
set(gca,'XScale','log','YScale','log','XLim',[.15 1.5],'YLim',[2 200])
cf = fit(vectorise(back_stim),cdata3.fA_gain,'power1','Upper',[Inf -1],'Lower',[-Inf -1]);
plot(back_stim,cf(back_stim),'r')
title('Gain estimation using pulses on background')
xlabel('Mean (background stimulus)')
ylabel('Gain (Hz/V)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Effect of stimulus offset
% What if we made an error of offset while measuring the stimulus? How does this affect our estimate of gain, and of Weber's Law? 

load('DA_model_fit_to_MSG.mat','p')
% also generate responses to simple steps 
cdata3 = cdata2;
back_stim = logspace(log10(min(mean_stim)),log10(max(mean_stim)),size(cdata.PID,2));
for i = 1:length(cdata.paradigm)
	S = ones(60e3,1)*back_stim(i);
	S(30e3:31e3) = 2*S(30e3:31e3);
	cdata3.PID(:,i) = S;
	R = DAModelv2(S,p(1));

	% compute gain
	dR = (max(R(30e3:31e3))-mean(R(28e3:29e3)));
	dS = back_stim(i);
	cdata3.fA_gain(i) = dR/dS;

end
figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(back_stim,cdata3.fA_gain,'k+')
set(gca,'XScale','log','YScale','log','XLim',[.15 1.5],'YLim',[2 200])
cf = fit(vectorise(back_stim),cdata3.fA_gain,'power1','Upper',[Inf -1],'Lower',[-Inf -1]);
plot(back_stim,cf(back_stim),'r')
title('Gain estimation using pulses on background')
xlabel('Mean (background stimulus)')
ylabel('Gain (Hz/V)')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% As we see from the figure, a 5% error in the offset of the stimulus can throw off our estimate of gain, and can introduce weird artifacts, including an apparent mismatch with Weber's Law (the DA model used here is the same, the only thing that differs is an offset in the stimulus). 

%% Effect of noise
% What is we add some noise to the data, and then attempt to extract fitlers? Here we introduce output noise in the responses and try to compute gain, and see how it affects our determination of Weber's Law.   

RandStream.setGlobalStream(RandStream('mt19937ar','Seed',1984)); 
noise = 100;

cdata2 = cdata;
for i = 1:width(cdata.PID)
	p(i).s0 = 0;
	cdata2.fA(:,i) = DAModelv2(cdata.PID(:,i),p(cdata.orn(i))) + filtfilt(ones(50,1),50,randn(60e3,1))*noise;
end
cdata2.LFP = NaN*cdata2.LFP;
cdata2 = cleanMSGdata(cdata2,'use_cache',false);

figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(mean_stim,cdata2.fA_gain,'k+')
set(gca,'XScale','log','YScale','log','XLim',[.15 1.5],'YLim',[2 200])
cf = fit(vectorise(mean_stim),cdata2.fA_gain,'power1','Upper',[Inf -1],'Lower',[-Inf -1]);
plot(mean_stim,cf(mean_stim),'r')
xlabel('Mean (white noise stimulus)')
ylabel('Gain (Hz/V)')
title('Gain estimation with additive noise')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% In this experiment, we added noise to all the trials, which led us to under-estimate the steepness of the gain dependency with mean stimulus. 

%% Metadata

pFooter


