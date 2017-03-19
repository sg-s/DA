

%% Fast gain control in naturalistic stimulus
% In this document, I attempt to use a model-free analysis method to demonstrate that there exists fast gain control in the naturalistic stimulus. 



pHeader;


%%
% In the following figure, I collect every whiff and measure the response of the neuron to that whiff. I then plot the response (LFP on the left, firing rate on the right), as a function of the whiff amplitude (x-axis) and the mean stimulus in the preceding 300ms (y-axis). The colours indicate the response magnitude (yellow = bigger). It's clear that as the whiff amplitude increases, the response amplitude increases. What is less clear is if the response amplitude decreases if the stimulus in the preceding 300 ms if large. 

cdata = consolidateData2(getPath(dataManager,'4608c42b12191b383c84fea52392ea97'));
[cdata, data] =  assembleScaledNatStim(cdata);

clear ws
for j = 1:3
	ws(j) = whiffStatistics(data(2).S(:,j),data(2).X(:,j),data(2).R(:,j),300,'MinPeakProminence',max(data(2).S(:,j)/100));
end

% plot the lfp
X = -vertcat(ws.peak_LFP);
R = vertcat(ws.peak_firing_rate);
S = vertcat(ws.stim_peaks);
Shat = vertcat(ws.stim_history_length_before_whiff);

rm_this = isnan(X) | isnan(R);
X(rm_this) = []; S(rm_this) = []; Shat(rm_this) = []; R(rm_this) = [];

% make a colors vector
c = parula(100);
cX = X - min(X); cX = cX/max(cX); cX = floor(1+cX*99);

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
scatter(S,Shat,256,cX,'filled');
set(gca,'XScale','log')
xlabel('Whiff amplitude (V)')
ylabel('\mu_{Stimulus} in preceding 300ms')
title('LFP responses')

% make a colors vector
c = parula(100);
cX = R - min(R); cX = cX/max(cX); cX = floor(1+cX*99);

subplot(1,2,2); hold on
scatter(S,Shat,256,cX,'filled');
set(gca,'XScale','log')
xlabel('Whiff amplitude (V)')
ylabel('\mu_{Stimulus} in preceding 300ms')
title('Firing rate responses')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end


%%
% To look at this more closely, I collect all responses to all whiffs, and then assign a color based on whether the response is bigger (yellow) or smaller (purple) than I would expect, by comparing it to similar-sized whiffs in the dataset.  In the following plot, I plot the responses, colour coded as indicated, as a function of the time since previous whiff, and the size of the previous whiff. LFP responses are on the left, and firing rate responses are on the right. The lines indicate the averages of the x & y values of all the dots along the x and y axes. Note that the mean of the purple cloud is always lower and to the right of the mean of the yellow cloud, meaning that when responses are smaller than we expect, the previous whiff was large, and occured earlier in the past, than when responses were bigger than we expect.  

S = data(2).S; 
X = data(2).X;
R = data(2).R;

clear ws
for j = 1:3
	ws(j) = whiffStatistics(data(2).S(:,j),data(2).X(:,j),data(2).R(:,j),300,'MinPeakProminence',max(data(2).S(:,j)/100));
end

for i = 1:2
	ws(i+1).stim_peak_loc = ws(i+1).stim_peak_loc + 70e3*i;
end

stim_peaks = vertcat(ws.stim_peaks);
stim_peak_loc = vertcat(ws.stim_peak_loc);
R_peak = vertcat(ws.peak_firing_rate);
X_peak = vertcat(ws.peak_LFP);

D_whiff_R  = NaN*stim_peaks; % the normalised deivations of the response to this whiff, from other
D_whiff_X  = NaN*stim_peaks;  
T_before = stim_peak_loc - circshift(stim_peak_loc,1); % the time to the preceding whiff
S_before = circshift(stim_peaks,1); % the peak of the preceding whiff


for i = 1:length(stim_peaks)
	set_responses_R = nonnans(R_peak(stim_peaks > stim_peaks(i)*(.9) & stim_peaks < stim_peaks(i)*(1.1)));
	set_responses_X = -nonnans(X_peak(stim_peaks > stim_peaks(i)*(.9) & stim_peaks < stim_peaks(i)*(1.1)));
	D_whiff_R(i) = R_peak(i)/mean(set_responses_R);
	D_whiff_X(i) = -X_peak(i)/mean(set_responses_X);
	
end

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on

% make a colors vector
c = parula(100);
cX = D_whiff_X - min(D_whiff_X); cX = cX/max(cX); cX = floor(1+cX*99);
cX(cX>50) = 100;
cX(cX<=50) = 1;

y = round(mean(T_before(cX>50)));
x = mean(S_before(cX>50));
plot([1e-4 1e4],[y y],'Color',c(100,:),'LineWidth',3)
plot([x x],[1 1e4],'Color',c(100,:),'LineWidth',3)

y = round(mean(T_before(cX<50)));
x = mean(S_before(cX<50));
plot([1e-4 1e4],[y y],'Color',c(1,:),'LineWidth',3)
plot([x x],[1 1e4],'Color',c(1,:),'LineWidth',3)


scatter(S_before(~isnan(cX)),T_before(~isnan(cX)),256,cX(~isnan(cX)),'filled');
xlabel('Previous whiff size (V)')
ylabel('Time since previous whiff (ms)')
title('Deviations in LFP (norm)')
set(gca,'YLim',[1e2 1e4],'YScale','log','XScale','log','XLim',[1e-2 10])


subplot(1,2,2); hold on

% make a colors vector
c = parula(100);
cX = D_whiff_R - min(D_whiff_R); cX = cX/max(cX); cX = floor(1+cX*99);
cX(cX>50) = 100;
cX(cX<100) = 1;

scatter(S_before,T_before,256,cX,'filled');

y = round(mean(T_before(cX>50)));
x = mean(S_before(cX>50));
plot([1e-4 1e4],[y y],'Color',c(100,:),'LineWidth',3)
plot([x x],[1 1e4],'Color',c(100,:),'LineWidth',3)

y = round(mean(T_before(cX<50)));
x = mean(S_before(cX<50));
plot([1e-4 1e4],[y y],'Color',c(1,:),'LineWidth',3)
plot([x x],[1 1e4],'Color',c(1,:),'LineWidth',3)


xlabel('Previous whiff size (V)')
ylabel('Time since previous whiff (ms)')
title('Deviations in firing rate (norm)')

set(gca,'YLim',[1e2 1e4],'YScale','log','XScale','log','XLim',[1e-2 10])

% canvas = axes;
% canvas.Position = [0 0 1 1];
% canvas.XTick = []; 
% canvas.YTick = [];
% uistack(canvas,'bottom');

% ch = colorbar(canvas,'east');
% caxis([min(D_whiff) max(D_whiff)]);
% ch.Position = [.92 .1 .01 .8];

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% To push this idea home, I plot the stimulus averaged over tau ms before the whiff for all whiffs that elicited a larger than expected response, vs. stimulus averaged over tau ms before the whiff for all whiffs that elicited a smaller than expected response. To be clear, each dot in this plot corresponds to a set of whiffs, all of similar (around 10%) size, that are split into two groups: one which elicited bigger than expected responses, and another that elicited smaller than expected responses. The mean of each of these groups defines a point in this plot. I then repeat this analysis for many different timescales.  

%%
% Note that most of the points are below the diagonal line at some timescales, suggesting that the stimulus before smaller than expected responses is LARGER than the stimulus before larger than expected responses. In other words, there is an inverse correlation between gain and the stimulus in the preceding window. 

S = data(2).S; 
X = data(2).X;
R = data(2).R;

clear ws
for j = 1:3
	ws(j) = whiffStatistics(data(2).S(:,j),data(2).X(:,j),data(2).R(:,j),300,'MinPeakProminence',max(data(2).S(:,j)/100));
end

for i = 1:2
	ws(i+1).stim_peak_loc = ws(i+1).stim_peak_loc + 70e3*i;
end
S = S(:);

stim_peaks = vertcat(ws.stim_peaks);
stim_peak_loc = vertcat(ws.stim_peak_loc);
R_peak = vertcat(ws.peak_firing_rate);
X_peak = vertcat(ws.peak_LFP);

window_size = [100 200 300 500 750 1e3 2e3 10e3];
buffer_size = 50;

clear plot_data
for oi = 1:length(window_size)

	small_X = NaN*stim_peaks;
	big_X = NaN*stim_peaks;

	small_R = NaN*stim_peaks;
	big_R = NaN*stim_peaks;

	for i = 1:length(stim_peaks)
		clear this_set
		this_set.R = (R_peak(stim_peaks > stim_peaks(i)*(.9) & stim_peaks < stim_peaks(i)*(1.1)));
		this_set.X = -(X_peak(stim_peaks > stim_peaks(i)*(.9) & stim_peaks < stim_peaks(i)*(1.1)));
		this_set.S_loc = (stim_peak_loc(stim_peaks > stim_peaks(i)*(.9) & stim_peaks < stim_peaks(i)*(1.1)));
		rm_this = isnan(this_set.R);
		this_set.R(rm_this) = [];
		this_set.X(rm_this) = [];
		this_set.S_loc(rm_this) = [];

		if length(this_set.R) > 3
			% cut out snippets from the stimulus 
			snippets = NaN(window_size(oi),length(this_set.R));
			for j = 1:length(this_set.R)
				try
					snippets(:,j) = S(this_set.S_loc(j)-window_size(oi)-buffer_size+1:this_set.S_loc(j)-buffer_size);
				catch
				end
			end
			big_R(i) = mean(mean(snippets(:,this_set.R > mean(this_set.R))));
			small_R(i) = mean(mean(snippets(:,this_set.R <= mean(this_set.R))));

			big_X(i) = mean(mean(snippets(:,this_set.X > mean(this_set.X))));
			small_X(i) = mean(mean(snippets(:,this_set.X <= mean(this_set.X))));
		end
	end
	plot_data(oi).big_R = big_R;
	plot_data(oi).small_R = small_R;

	plot_data(oi).big_X = big_X;
	plot_data(oi).small_X = small_X;
end

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:length(plot_data)
	autoPlot(length(plot_data),i); hold on
	plot(plot_data(i).small_R,plot_data(i).big_R,'k.','MarkerSize',20)
	set(gca,'XLim',[1e-3 1],'YLim',[1e-3 1],'XScale','log','YScale','log','XTick',[1e-3 1e-2 1e-1 1e0])
	plot([1e-6 0.7],[1e-6 0.7],':','Color',[.5 .5 .5])
	title(['\tau = ' oval(window_size(i)) ' ms'])
	if i == 1
		xlabel(['Stimulus \tau ms before smaller' char(10) 'than expected spiking'])
		ylabel(['Stimulus \tau ms before larger' char(10) 'than expected spiking'])
	end
end

prettyFig();

if being_published
	snapnow
	delete(gcf)
end



%%
% Now I do the same thing for the LFP. 

figure('outerposition',[0 0 1400 800],'PaperUnits','points','PaperSize',[1400 800]); hold on
for i = 1:length(plot_data)
	autoPlot(length(plot_data),i); hold on
	plot(plot_data(i).small_X,plot_data(i).big_X,'k.','MarkerSize',20)
	set(gca,'XLim',[1e-3 1],'YLim',[1e-3 1],'XScale','log','YScale','log','XTick',[1e-3 1e-2 1e-1 1e0])
	plot([1e-6 0.7],[1e-6 0.7],':','Color',[.5 .5 .5])
	title(['\tau = ' oval(window_size(i)) ' ms'])

	if i == 1
		xlabel(['Stimulus \tau ms before smaller' char(10) 'than expected LFP'])
		ylabel(['Stimulus \tau ms before larger' char(10) 'than expected LFP'])
	end
end


prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Now I plot the ratio of the stimulus before bigger than expected and smaller than expected responses as a function of the timescale we look over. This gives us a sense of how strong this effect is, and over what timescale it operates over. Note that the precise form of this plot is heavily dependent on the stimulus that we use. 

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
x = window_size;
y = NaN*window_size;
ye = NaN*window_size;
for i = 1:length(window_size)
	temp = nonnans(plot_data(i).small_X./plot_data(i).big_X);
	y(i) = mean(temp);
	ye(i) = sem(temp); 
end
errorbar(x,y,ye)
xlabel('\tau (ms)')
ylabel(['Stimulus before smaller LFP/' char(10) 'stimulus before bigger LFP'])
set(gca,'XScale','log','XTick',[10 100 1e3 1e4 1e5],'XLim',[100 1.1e4])

subplot(1,2,2); hold on
x = window_size;
y = NaN*window_size;
ye = NaN*window_size;
for i = 1:length(window_size)
	temp = nonnans(plot_data(i).small_R./plot_data(i).big_R);
	y(i) = mean(temp);
	ye(i) = sem(temp); 
end
errorbar(x,y,ye)
xlabel('\tau (ms)')
ylabel(['Stimulus before smaller spiking/' char(10) 'stimulus before bigger spiking'])
set(gca,'XScale','log','XTick',[10 100 1e3 1e4 1e5],'XLim',[100 1.1e4])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% Is there a way to determine what in the preceding stimulus causes these responses to vary from the mean? To do this, I collect all whiffs, and try to back out a linear predictor that uses the past stimulus to determine whether the responses to this whiff will be bigger or smaller than expected. In the following figure, I back out these filters, for the LFP and the firing rate. From the poor goodness of fits, it doesn't look like this works. 

stim_peaks = vertcat(ws.stim_peaks);
stim_peak_loc = vertcat(ws.stim_peak_loc);
R_peak = vertcat(ws.peak_firing_rate);
X_peak = vertcat(ws.peak_LFP);
X_D = NaN*X_peak;
R_D = NaN*X_peak;
Shat = NaN(10,length(X_peak)); % 10 100ms chunks in the past
S = S - min(S);

% compute deviations from expected values for each whiff
for i = 1:length(stim_peaks)
	clear this_set
	idx = stim_peaks > stim_peaks(i)*(.9) & stim_peaks < stim_peaks(i)*(1.1);
	this_set.R = (R_peak(idx));
	this_set.X = -(X_peak(idx));
	this_set.S_loc = (stim_peak_loc(idx));

	rm_this = isnan(this_set.R);
	this_set.R(rm_this) = [];
	this_set.X(rm_this) = [];
	this_set.S_loc(rm_this) = [];

	if length(this_set.R) > 3
		R_D(i) = R_peak(i)/mean(this_set.R) - 1; 
		X_D(i) = X_peak(i)/mean(this_set.X) - 1; 
		% also assemble the Shat vector
		if stim_peak_loc(i) > size(Shat,1)*100
			S_snippet = S(stim_peak_loc(i)-size(Shat,1)*100+1:stim_peak_loc(i));
			Shat(:,i) = interp1(1:length(S_snippet),S_snippet,linspace(1,length(S_snippet),size(Shat,1)));
		end
	end
end

% remove whiffs for which we could not estimate the deviation because we didn't have enough similarly sized whiffs 
rm_this = isnan(sum(Shat)) | isnan(R_D)';
S_ = Shat(:,~rm_this);
C = S_*S_';
C = C + eye(length(C))*mean(eig(C));

figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(1,2,1); hold on
D_ = X_D(~rm_this);
K = (C\(S_*D_));
bar((1:length(K))/10,flipud(K))
xlabel('Lag (s)')
ylabel('LFP filter')
title(['r^2 = ' oval(rsquare(S_'*K,D_))])

subplot(1,2,2); hold on
D_ = R_D(~rm_this);
K = (C\(S_*D_));
bar((1:length(K))/10,flipud(K))
xlabel('Lag (s)')
ylabel('Firing filter')
title(['r^2 = ' oval(rsquare(S_'*K,D_))])

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%%
% This is annoying: while the LFP filter looks "nice", in that it is all negative, the r^2 is 0, meaning this filter has no predictive power. And while the firing rate filter has a nice minimum at ~300ms, it is hard to rule out the contribution of a response filter with a negative lobe from this. 

%% Analysis pipeline 
% Now I describe the analysis pipeline I developed to investigate whether variations in response arise from variations in the pulse height, or from the history of the stimulus. This method of analysis is so fiendishly complicated that I made a handy graphic that explains what I'm doing: 

figure('outerposition',[0 0 1000 655],'PaperUnits','points','PaperSize',[1000 655]); hold on
% show the cartoons
o = imread('../images/fast_gain_control_analysis.png');
imagesc(o);
axis ij
axis image
axis off

prettyFig();

if being_published
	snapnow
	delete(gcf)
end




;;;;;;;;  ;;;;;;;;    ;;;    ;;          ;;;;;;;;     ;;;    ;;;;;;;;    ;;;    
;;     ;; ;;         ;; ;;   ;;          ;;     ;;   ;; ;;      ;;      ;; ;;   
;;     ;; ;;        ;;   ;;  ;;          ;;     ;;  ;;   ;;     ;;     ;;   ;;  
;;;;;;;;  ;;;;;;   ;;     ;; ;;          ;;     ;; ;;     ;;    ;;    ;;     ;; 
;;   ;;   ;;       ;;;;;;;;; ;;          ;;     ;; ;;;;;;;;;    ;;    ;;;;;;;;; 
;;    ;;  ;;       ;;     ;; ;;          ;;     ;; ;;     ;;    ;;    ;;     ;; 
;;     ;; ;;;;;;;; ;;     ;; ;;;;;;;;    ;;;;;;;;  ;;     ;;    ;;    ;;     ;; 

%% Real data: ab2A responses to 2-butanone
% In this section I apply this analysis method to ab2A responses to naturalistic stimuli constructed with 2-butanone. Since this stimulus is presented many times in scaled versions of itself, it covers the entire dynamic range of the neuron. (a-c) Analysis of fast gain control in LFP. (a) LFP responses vs. whiff amplitude, for sets of whiffs with similar amplitude (see methods in previous figure). The variation in whiff amplitude cannot describe the variation in the LFP responses. Lines connect whiffs of the same size that come from the same set. (b) LFP responses vs. mean stimulus in the preceding 300ms. The strong negative correlation shows what we saw in the raw data: that LFP responses to whiffs following other whiffs are smaller than LFp responses to whiffs that occur in isolation. (c) Spearman correlation between LFP responses and whiff amplitude (black) and between LFP responses and mean stimulus in the preceding window. Crosses indicate significance (p<0.05). Note that there is a pronounced dip at around ~300ms. 

%%
% (d-f) Similar analysis, but for the firing rates. 

[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,data(2).X,data(2).R);

figure('outerposition',[0 0 1400 901],'PaperUnits','points','PaperSize',[1400 901]); hold on
subplot(2,3,1); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,x))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('LFP response (norm)')

subplot(2,3,2); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,x))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('LFP response (norm)')

subplot(2,3,4); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,r))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('Firing response (norm)')

subplot(2,3,5); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,r))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('Firing response (norm)')

% now vary the stimulus history length 
all_history_lengths = unique(floor(logspace(1,4,60)));

rho_x_whiff = NaN*all_history_lengths;
rho_x_gain = NaN*all_history_lengths;
rho_r_whiff = NaN*all_history_lengths;
rho_r_gain = NaN*all_history_lengths;

rho_x_whiff_p = NaN*all_history_lengths;
rho_x_gain_p = NaN*all_history_lengths;
rho_r_whiff_p = NaN*all_history_lengths;
rho_r_gain_p = NaN*all_history_lengths;

for i = 1:length(all_history_lengths)
	[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,data(2).X,data(2).R,'history_length',all_history_lengths(i));

	% compute correlations 
	[rho_x_gain(i), rho_x_gain_p(i)] = corr(s,x,'type','Spearman');
	[rho_x_whiff(i), rho_x_whiff_p(i)] = corr(whiff_s,x,'type','Spearman');
	[rho_r_gain(i), rho_r_gain_p(i)] = corr(s,r,'type','Spearman');
	[rho_r_whiff(i), rho_r_whiff_p(i)] = corr(whiff_s,r,'type','Spearman');


end

% emphasise significant points
subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k');
plot(all_history_lengths,rho_x_gain,'r');


subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k');
plot(all_history_lengths,rho_r_gain,'r');



rho_x_gain(rho_x_gain_p>.05) = NaN;
rho_x_whiff(rho_x_whiff_p>.05) = NaN;
rho_r_gain(rho_r_gain_p>.05) = NaN;
rho_r_whiff(rho_r_whiff_p>.05) = NaN;

subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k+');
plot(all_history_lengths,rho_x_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')

subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k+');
plot(all_history_lengths,rho_r_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')


prettyFig();
labelFigure('x_offset',0)

if being_published
	snapnow
	delete(gcf)
end




;;     ;;    ;;;    ;;       ;;;; ;;;;;;;;     ;;;    ;;;;;;;; ;;;;  ;;;;;;;  ;;    ;; 
;;     ;;   ;; ;;   ;;        ;;  ;;     ;;   ;; ;;      ;;     ;;  ;;     ;; ;;;   ;; 
;;     ;;  ;;   ;;  ;;        ;;  ;;     ;;  ;;   ;;     ;;     ;;  ;;     ;; ;;;;  ;; 
;;     ;; ;;     ;; ;;        ;;  ;;     ;; ;;     ;;    ;;     ;;  ;;     ;; ;; ;; ;; 
 ;;   ;;  ;;;;;;;;; ;;        ;;  ;;     ;; ;;;;;;;;;    ;;     ;;  ;;     ;; ;;  ;;;; 
  ;; ;;   ;;     ;; ;;        ;;  ;;     ;; ;;     ;;    ;;     ;;  ;;     ;; ;;   ;;; 
   ;;;    ;;     ;; ;;;;;;;; ;;;; ;;;;;;;;  ;;     ;;    ;;    ;;;;  ;;;;;;;  ;;    ;; 

;;    ;; ;;          ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;       
;;;   ;; ;;          ;;;   ;;; ;;     ;; ;;     ;; ;;       ;;       
;;;;  ;; ;;          ;;;; ;;;; ;;     ;; ;;     ;; ;;       ;;       
;; ;; ;; ;;          ;; ;;; ;; ;;     ;; ;;     ;; ;;;;;;   ;;       
;;  ;;;; ;;          ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;   ;;; ;;          ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;    ;; ;;;;;;;;    ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;;;;;;; 


%% Validation: NL model (LFP)
% As a validation, I repeat this analysis using synthetic data generated by a NL model. (a-c) Analysis using synthetic data generated by a NL model fit to the LFP responses. (d-f) Analysis using real LFP responses. (a) NL model responses vs. whiff amplitude, for similar sized whiffs. Whiff amplitudes cannot account for observed variation in responses. (b) NL model responses vs. mean stimulus in preceding 300ms. The mean stimulus in the preceding 300ms can account for the variation in the synthetic data, but it goes the "wrong way" -- NL model responses tend to be bigger when the mean stimulus in the preceding 300ms is larger -- consistent with an integrating, but non-adapting system. (c) Spearman correlation for NL model responses vs. whiff amplitude (black) and vs. mean stimulus in preceding window (red). Crossed indicate significant correlation. Note that for this synthetic data set, the mean stimulus in any preceding window always tends to amplify responses. 

%%
% (e-f) as in the previous figure. 

clear p
p.k_D = 0.21246;
p.n = 1.4688;


% generate responses using this model 
warning off
i = 2;
for j = 1:size(data(i).X,2)
	data(i).P(:,j) = NLModel([data(i).S(:,j) - min(data(i).S(:,j)) data(i).X(:,j)] ,p);
end

warning on

% we now attempt to estimate a gain control timescale from the data
% first collect all whiffs
[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,data(2).P,-data(2).X);

figure('outerposition',[0 0 1400 901],'PaperUnits','points','PaperSize',[1400 901]); hold on
subplot(2,3,1); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,x))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('NL model response (norm)')

subplot(2,3,2); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,x))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('NL model response (norm)')

subplot(2,3,4); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,r))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('LFP response (norm)')

subplot(2,3,5); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,r))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('LFP response (norm)')

% now vary the stimulus history length 
all_history_lengths = unique(floor(logspace(1,4,60)));

rho_x_whiff = NaN*all_history_lengths;
rho_x_gain = NaN*all_history_lengths;
rho_r_whiff = NaN*all_history_lengths;
rho_r_gain = NaN*all_history_lengths;

rho_x_whiff_p = NaN*all_history_lengths;
rho_x_gain_p = NaN*all_history_lengths;
rho_r_whiff_p = NaN*all_history_lengths;
rho_r_gain_p = NaN*all_history_lengths;

for i = 1:length(all_history_lengths)
	[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,data(2).P,-data(2).X,'history_length',all_history_lengths(i));

	% compute correlations 
	[rho_x_gain(i), rho_x_gain_p(i)] = corr(s,x,'type','Spearman');
	[rho_x_whiff(i), rho_x_whiff_p(i)] = corr(whiff_s,x,'type','Spearman');
	[rho_r_gain(i), rho_r_gain_p(i)] = corr(s,r,'type','Spearman');
	[rho_r_whiff(i), rho_r_whiff_p(i)] = corr(whiff_s,r,'type','Spearman');


end

% emphasise significant points
subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k');
plot(all_history_lengths,rho_x_gain,'r');


subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k');
plot(all_history_lengths,rho_r_gain,'r');



rho_x_gain(rho_x_gain_p>.05) = NaN;
rho_x_whiff(rho_x_whiff_p>.05) = NaN;
rho_r_gain(rho_r_gain_p>.05) = NaN;
rho_r_whiff(rho_r_whiff_p>.05) = NaN;

subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k+');
plot(all_history_lengths,rho_x_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')

subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k+');
plot(all_history_lengths,rho_r_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')


prettyFig();
labelFigure('x_offset',0);

if being_published
	snapnow
	delete(gcf)
end


;;     ;;    ;;;    ;;       ;;;; ;;;;;;;;     ;;;    ;;;;;;;; ;;;;  ;;;;;;;  ;;    ;; 
;;     ;;   ;; ;;   ;;        ;;  ;;     ;;   ;; ;;      ;;     ;;  ;;     ;; ;;;   ;; 
;;     ;;  ;;   ;;  ;;        ;;  ;;     ;;  ;;   ;;     ;;     ;;  ;;     ;; ;;;;  ;; 
;;     ;; ;;     ;; ;;        ;;  ;;     ;; ;;     ;;    ;;     ;;  ;;     ;; ;; ;; ;; 
 ;;   ;;  ;;;;;;;;; ;;        ;;  ;;     ;; ;;;;;;;;;    ;;     ;;  ;;     ;; ;;  ;;;; 
  ;; ;;   ;;     ;; ;;        ;;  ;;     ;; ;;     ;;    ;;     ;;  ;;     ;; ;;   ;;; 
   ;;;    ;;     ;; ;;;;;;;; ;;;; ;;;;;;;;  ;;     ;;    ;;    ;;;;  ;;;;;;;  ;;    ;; 

;;    ;; ;;       ;;    ;;    ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;       
;;;   ;; ;;       ;;;   ;;    ;;;   ;;; ;;     ;; ;;     ;; ;;       ;;       
;;;;  ;; ;;       ;;;;  ;;    ;;;; ;;;; ;;     ;; ;;     ;; ;;       ;;       
;; ;; ;; ;;       ;; ;; ;;    ;; ;;; ;; ;;     ;; ;;     ;; ;;;;;;   ;;       
;;  ;;;; ;;       ;;  ;;;;    ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;   ;;; ;;       ;;   ;;;    ;;     ;; ;;     ;; ;;     ;; ;;       ;;       
;;    ;; ;;;;;;;; ;;    ;;    ;;     ;;  ;;;;;;;  ;;;;;;;;  ;;;;;;;; ;;;;;;;; 


%% Validation: NLN model (firing rate)
% As a validation, I repeat this analysis using synthetic data generated by a NLN model. This NLN model has a negative lobe, since it is fit to the firing rate data, and this could in principle give us something that looks like "fast gain control".  (a-c) Analysis using synthetic data generated by a NLN model fit to the firing rate responses. (d-f) Analysis using real firing rate responses. (a) NLN model responses vs. whiff amplitude, for similar sized whiffs. Whiff amplitudes cannot account for observed variation in responses. (b) NLN model responses vs. mean stimulus in preceding 300ms. The mean stimulus in the preceding 300ms can account for the variation in the synthetic data, similar to what we see in the data. (c) Spearman correlation for NLN model responses vs. whiff amplitude (black) and vs. mean stimulus in preceding window (red). Crossed indicate significant correlation. Note that for this synthetic data set, the mean stimulus in any preceding window always tends to amplify responses. 

%%
% (d-f) as in the previous figure with real data. 

clear p
p.  k_D = 0.1109;
p.    n = .9688;

% generate responses using this model 
warning off
i = 2;
for j = 1:size(data(i).X,2)
	try
		data(i).P(:,j) = NLNmodel([data(i).S(:,j) - min(data(i).S(:,j)) data(i).R(:,j)] ,p);
	catch
	end
end

warning on

% we now attempt to estimate a gain control timescale from the data
% first collect all whiffs
[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,-data(2).P,data(2).R);

figure('outerposition',[0 0 1400 901],'PaperUnits','points','PaperSize',[1400 901]); hold on
subplot(2,3,1); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,x))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('NLN model response (norm)')

subplot(2,3,2); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = x(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,x))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('NLN model response (norm)')

subplot(2,3,4); hold on
for i = 1:max(group)
	tempx = whiff_s(group==i)/mean(whiff_s(group==i));
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(whiff_s,r))],'Location','northeast')
xlabel('Whiff amplitude (norm)')
ylabel('Firing response (norm)')

subplot(2,3,5); hold on
for i = 1:max(group)
	tempx = s(group==i);
	tempy = r(group==i);
	[tempx,idx] = sort(tempx);
	plot(tempx,tempy(idx),'+-');
end
l = plot(NaN,NaN,'k+');
legend(l,['\rho = ' oval(spear(s,r))],'Location','northeast')
xlabel('\mu_{Stimulus} in preceding 300ms')
ylabel('Firing response (norm)')

% now vary the stimulus history length 
all_history_lengths = unique(floor(logspace(1,4,60)));

rho_x_whiff = NaN*all_history_lengths;
rho_x_gain = NaN*all_history_lengths;
rho_r_whiff = NaN*all_history_lengths;
rho_r_gain = NaN*all_history_lengths;

rho_x_whiff_p = NaN*all_history_lengths;
rho_x_gain_p = NaN*all_history_lengths;
rho_r_whiff_p = NaN*all_history_lengths;
rho_r_gain_p = NaN*all_history_lengths;

for i = 1:length(all_history_lengths)
	[s,x,r,whiff_s,group] = binnedWhiffGainAnalysis(data(2).S,-data(2).P,data(2).R,'history_length',all_history_lengths(i));

	% compute correlations 
	[rho_x_gain(i), rho_x_gain_p(i)] = corr(s,x,'type','Spearman');
	[rho_x_whiff(i), rho_x_whiff_p(i)] = corr(whiff_s,x,'type','Spearman');
	[rho_r_gain(i), rho_r_gain_p(i)] = corr(s,r,'type','Spearman');
	[rho_r_whiff(i), rho_r_whiff_p(i)] = corr(whiff_s,r,'type','Spearman');


end

% emphasise significant points
subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k');
plot(all_history_lengths,rho_x_gain,'r');


subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k');
plot(all_history_lengths,rho_r_gain,'r');



rho_x_gain(rho_x_gain_p>.05) = NaN;
rho_x_whiff(rho_x_whiff_p>.05) = NaN;
rho_r_gain(rho_r_gain_p>.05) = NaN;
rho_r_whiff(rho_r_whiff_p>.05) = NaN;

subplot(2,3,3); hold on
plot(all_history_lengths,rho_x_whiff,'k+');
plot(all_history_lengths,rho_x_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')

subplot(2,3,6); hold on
plot(all_history_lengths,rho_r_whiff,'k+');
plot(all_history_lengths,rho_r_gain,'r+');
set(gca,'XScale','log')
xlabel('Stimulus history length (ms)')
ylabel('\rho')


prettyFig();

labelFigure('x_offset',0)

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;


