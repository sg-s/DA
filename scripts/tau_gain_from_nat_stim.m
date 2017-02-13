


pHeader;

%% Estimating the timescale of gain control from naturalistic stimuli
% In this document, I attempt to estimate the gain control timescale from naturalistic stimuli, where brief whiff of odorant are presented to neurons, with widely varying amplitudes, some in isolation, and some in dense whiffs where other whiffs precede them. 

%% Background
% I attempted to this previously by looking at how deviations from a linear model varied with the mean stimulus in some preceding window. However, there was a fatal flaw in this analysis: these deviations could be largely accounted for by an input nonlinearity, that correlations in the stimulus could then manifest as "dynamic" gain control. 

%%
% In this document, I work around this problem by first fitting, explicitly, a nonlinear-linear model to the data. These models reproduce the data very well (with a very high $r^2$). Now that I have accounted for the input nonlinearity, I study how deviations of the neuron response from this NLN model correlate with the mean stimulus in some preceding window. 

   ;;;    ;;;;;;;;   ;;;;;;;     ;;;    
  ;; ;;   ;;     ;; ;;     ;;   ;; ;;   
 ;;   ;;  ;;     ;;        ;;  ;;   ;;  
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;; 
;;;;;;;;; ;;     ;; ;;        ;;;;;;;;; 
;;     ;; ;;     ;; ;;        ;;     ;; 
;;     ;; ;;;;;;;;  ;;;;;;;;; ;;     ;; 

 ;;;;;;;          ;;;;;;;;  ;;     ;; ;;;;;;;;    ;;;    ;;    ;;  ;;;;;;;  ;;    ;; ;;;;;;;; 
;;     ;;         ;;     ;; ;;     ;;    ;;      ;; ;;   ;;;   ;; ;;     ;; ;;;   ;; ;;       
       ;;         ;;     ;; ;;     ;;    ;;     ;;   ;;  ;;;;  ;; ;;     ;; ;;;;  ;; ;;       
 ;;;;;;;  ;;;;;;; ;;;;;;;;  ;;     ;;    ;;    ;;     ;; ;; ;; ;; ;;     ;; ;; ;; ;; ;;;;;;   
;;                ;;     ;; ;;     ;;    ;;    ;;;;;;;;; ;;  ;;;; ;;     ;; ;;  ;;;; ;;       
;;                ;;     ;; ;;     ;;    ;;    ;;     ;; ;;   ;;; ;;     ;; ;;   ;;; ;;       
;;;;;;;;;         ;;;;;;;;   ;;;;;;;     ;;    ;;     ;; ;;    ;;  ;;;;;;;  ;;    ;; ;;;;;;;; 


%% ab2A and 2-butanone
% The first dataset I will try this on is with ab2A and 2-butanone. The data we have here is very broadly distributed, as we have the same naturalistic stimulus at three different scales. 



%% 
% Can we see if there is any gain control in this data? In the following figure, I first fit NLN-models neuron-by-neuron to the data. The NLN models have only two parameters (the $k_D$ and the $n$) that are fit parametrically. Other parameters (like the filter) are fit non-parametrically. In the following figure, I plot the results of each neuron in a separate colour. The first plot shows the distribution of the deviations of the NLN model predictions from the measured response. The second plot shows the Spearman correlation between the deviations and the mean stimulus in some preceding window, as a function of window length. Note that these plots tend to have a minimum at some defined timescale. The dotted and dashed lines indicate the autocorrelation times for the stimulus and the response respectively. 


% get all data 
cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[~, data] =  assembleScaledNatStim(cdata);

this_orn = 2;
clear fd
for i = 1:size(data(this_orn).X,2)
	S = data(this_orn).S(:,i); S = S - min(S);
	R = -data(this_orn).X(:,i);
	fd(i).stimulus = [S R];
	fd(i).response = R;
end
fd(any(sum(isnan(horzcat(fd.response))))) = [];

clear p

p(1).k_D = .6578;
p(1).n = .7188;

p(2).  k_D = 0.1109;
p(2).    n = .9688;

p(3).k_D = .1187;
p(3).n = .7812;

% generate responses using this model 
warning off
for i = 1:length(data)
	for j = 1:size(data(i).X,2)
		try
			data(i).P(:,j) = NLNmodel([data(i).S(:,j) - min(data(i).S(:,j)) data(i).R(:,j)] ,p(i));
		catch
		end
	end
end
warning on

figure('outerposition',[0 0 1220 601],'PaperUnits','points','PaperSize',[1220 601]); hold on
clear ax
ax(2) = subplot(1,2,1); hold on
ax(4) = subplot(1,2,2); hold on

c = lines(3);

for i = 1:3
	plot_tau_gain_nat_stim(data(i),ax,c(i,:));
end

prettyFig();
suptitle('ab2A -- 2-butanone')

if being_published
	snapnow
	delete(gcf)
end



   ;;;    ;;;;;;;;   ;;;;;;;     ;;;    
  ;; ;;   ;;     ;; ;;     ;;   ;; ;;   
 ;;   ;;  ;;     ;;        ;;  ;;   ;;  
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;; 
;;;;;;;;; ;;     ;;        ;; ;;;;;;;;; 
;;     ;; ;;     ;; ;;     ;; ;;     ;; 
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;; 


 ;;;;;;;          ;;;;;;;;  ;;     ;; ;;;;;;;;    ;;;    ;;    ;;  ;;;;;;;  ;;    ;; ;;;;;;;; 
;;     ;;         ;;     ;; ;;     ;;    ;;      ;; ;;   ;;;   ;; ;;     ;; ;;;   ;; ;;       
       ;;         ;;     ;; ;;     ;;    ;;     ;;   ;;  ;;;;  ;; ;;     ;; ;;;;  ;; ;;       
 ;;;;;;;  ;;;;;;; ;;;;;;;;  ;;     ;;    ;;    ;;     ;; ;; ;; ;; ;;     ;; ;; ;; ;; ;;;;;;   
;;                ;;     ;; ;;     ;;    ;;    ;;;;;;;;; ;;  ;;;; ;;     ;; ;;  ;;;; ;;       
;;                ;;     ;; ;;     ;;    ;;    ;;     ;; ;;   ;;; ;;     ;; ;;   ;;; ;;       
;;;;;;;;;         ;;;;;;;;   ;;;;;;;     ;;    ;;     ;; ;;    ;;  ;;;;;;;  ;;    ;; ;;;;;;;; 


%% ab3A and 2-butanone.
% I repeat the same analysis for ab3A and 2-butanone. Since I only have only neuron here, I also plot some extra things like a comparison of the response to the prediction, and a plot of the deviations vs. the mean stimulus in the preceding window.  

% get all data 
cdata = consolidateData2(getPath(dataManager,'c2bce18a6b0a7e89e9c6832dcc27e39b'));
[~, data] =  assembleScaledNatStim(cdata);



clear p
 p.k_D = 3.5894;
 p.  n = 0.5625;


% generate responses using this model 
clear K
for j = 1:size(data.X,2)
	data.P(:,j) = NLNmodel([data.S(:,j) - min(data.S(:,j)) data.R(:,j)] ,p);
end

figure('outerposition',[0 0 1000 901],'PaperUnits','points','PaperSize',[1000 901]); hold on
clear ax
for i = 1:4
	ax(i) = subplot(2,2,i); hold on
end

plot_tau_gain_nat_stim(data,ax);

prettyFig();
suptitle('ab3A -- 2-butanone')

if being_published
	snapnow
	delete(gcf)
end




   ;;;    ;;;;;;;;   ;;;;;;;     ;;;        ;;;;;;;     ;;;     ;;;;;;  
  ;; ;;   ;;     ;; ;;     ;;   ;; ;;      ;;     ;;   ;; ;;   ;;    ;; 
 ;;   ;;  ;;     ;;        ;;  ;;   ;;            ;;  ;;   ;;  ;;       
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;;     ;;;;;;;  ;;     ;; ;;       
;;;;;;;;; ;;     ;;        ;; ;;;;;;;;;    ;;        ;;;;;;;;; ;;       
;;     ;; ;;     ;; ;;     ;; ;;     ;;    ;;        ;;     ;; ;;    ;; 
;;     ;; ;;;;;;;;   ;;;;;;;  ;;     ;;    ;;;;;;;;; ;;     ;;  ;;;;;;  

%% ab3A and ethyl acetate
% In this section, I repeat the analysis but with ethyl acetate naturalistic stimulus applied to ab3A neurons. 

% load the data
cdata = consolidateKontrollerData(getPath(dataManager,'f70e37a7db469b88c0fc79ff5e828e9d'));
v2struct(cdata); clear cdata

clear data
for i = 1:max(orn)
	R = fA(:,orn==i);
	S = PID(:,orn==i);	
	rm_this = sum(R)==0;
	data(i).S = nanmean(S(:,~rm_this),2);
	data(i).S = data(i).S - min(data(i).S);
	data(i).R = nanmean(R(:,~rm_this),2);
end

% fit a NLN model to this
% for i = 1:length(data)
% 	fd(i).stimulus = [data(i).S data(i).R];
% 	fd(i).response = data(i).R;
% end

% for i = 1:length(fd)
% 	p(i) = fitModel2Data(@NLNmodel,fd(i),'p0',p(i),'nsteps',10);
% end

load('.cache/NLN_model_fit_to_ab3A_2ac_nat_stim.mat','p')

% use these parameters to make NLN model predictions
for i = 1:max(orn)
	data(i).P = NLNmodel([data(i).S data(i).R] ,p(i));
end

figure('outerposition',[0 0 1201 600],'PaperUnits','points','PaperSize',[1201 600]); hold on
clear ax
ax(4) = subplot(1,2,2); hold on
ax(2) = subplot(1,2,1); hold on

c = lines(length(data));

for i = [1 2 4 5 6]
	plot_tau_gain_nat_stim(data(i),ax,c(i,:));
end

prettyFig();
suptitle('ab3A -- ethyl acetate')

if being_published
	snapnow
	delete(gcf)
end



%% Version Info
%
pFooter;


