% 
% 
% created by Srinivas Gorur-Shandilya at 4:36 , 21 March 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


pHeader;




%% Idealised Spike Filters
% In this document, I explore the possibility that the true shape of the LFP to firing rate filter is narrow, and causal as opposed to broad and a-causal, as we get from our raw data. The reasoning is that since the LFP signal is correlated, and not necessarily Gaussian, our filter estimation techniques aren't guaranteed to work. This follows from Rachel Wilson's paper where they show a idealised, causal, differentiating filter does almost as good of a job as their data-constrained filter.  

% get the data
[PID, LFP, fA, paradigm,~, ~, AllControlParadigms] = consolidateData('/local-data/DA-paper/LFP-MSG/september',1);


% sort the paradigms sensibly
sort_value = [];
for i = 1:length(AllControlParadigms)
	sort_value(i) = (mean(AllControlParadigms(i).Outputs(1,:)));
end
[~,idx] = sort(sort_value);

AllControlParadigms = AllControlParadigms(idx);
paradigm_new = paradigm*NaN;
for i = 1:length(idx)
	paradigm_new(paradigm == idx(i)) = i;
end
paradigm = paradigm_new;

% remove baseline from all PIDs
for i = 1:width(PID)
	PID(:,i) = PID(:,i) - mean(PID(1:5e3,i));
end

% throw out trials where we didn't record the LFP, for whatever reason
not_LFP = find((max(abs(LFP))) < 0.1);
LFP(:,not_LFP) = NaN;

% throw our bad traces
bad_trials = (sum(fA) == 0 | isnan(sum(fA)) |  isnan(sum(LFP)));
LFP(:,bad_trials) = [];
PID(:,bad_trials) = [];
fA(:,bad_trials) = [];
paradigm(bad_trials) = [];

% band pass all the LFP
try 
	load('/local-data/DA-paper/LFP-MSG/september/filtered_LFP.mat','filtered_LFP')
catch
	filtered_LFP = LFP;
	for i = 1:width(LFP)
		filtered_LFP(:,i) = 10*bandPass(LFP(:,i),1e4,Inf);
	end
end

% define limits on data
a = 10e3; z = 50e3;

% extract filters and compute gains
[K1,K1p,K1_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);
[K2,K2p,K2_gain] = extractFilters(filtered_LFP,fA,'use_cache',true,'a',a,'z',z);
[K3,K3p,K3_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);
ft = 1e-3*(1:length(K1)) - .1;

%%
% In the following figure, we plot the filters backed out using standard methods together with the idealised filter that is force to be causal. We also plot the predictions using each of these filters and compare how well they predict the response. Each row shows a different trial from the full dataset. 

figure('outerposition',[0 0 800 800],'PaperUnits','points','PaperSize',[800 800]); hold on
do_these = find(paradigm==1);
for i = 1:3
	subplot(3,3,3*(i-1)+1), hold on
	d.stimulus = filtered_LFP(a:z,do_these(i));
	d.response = fA(a:z,do_these(i));
	p = fitModel2Data(@idealisedFilterModel,d,'nsteps',0);
	[R,K] = idealisedFilterModel(filtered_LFP(:,do_these(i)),p);
	plot(ft,K2(:,do_these(i))/norm(K2(:,do_these(i))),'k'), hold on
	plot(1e-3*(1:1e3),K/norm(K),'r')
	set(gca,'XLim',[-.2 .6])


	subplot(3,3,3*(i-1)+2), hold on
	plot(K2p(a:10:z,do_these(i)),fA(a:10:z,do_these(i)),'k.')
	legend(['r^2 = ' oval(rsquare(K2p(a:z,do_these(i)),fA(a:z,do_these(i))))],'Location','southeast')
	xlabel('Projection using normal filter')
	ylabel('Firing Rate (Hz)')

	subplot(3,3,3*(i-1)+3), hold on
	plot(R(a:10:z),fA(a:10:z,do_these(i)),'r.')
	legend(['r^2 = ' oval(rsquare(R(a:z),fA(a:z,do_these(i))))],'Location','southeast')
	xlabel('Projection using idealised filter')
	ylabel('Firing Rate (Hz)')
end
prettyFig('fs',12)

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;

