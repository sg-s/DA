%% fig_LFP.m
% makes figure showing how gain changes at the LFP level, and from the LFP to the firing machinery
% 


pHeader;

% uses dataManager for data integrity
dm = dataManager;

opacity = .5;

figure('outerposition',[0 0 1600 720],'PaperUnits','points','PaperSize',[1600 720]); hold on


[PID, LFP, fA, paradigm,~, ~, AllControlParadigms] = consolidateData(dm.getPath('bf79dfd769a97089e42beb0660174e84'),1);
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
	load(dm.getPath('213e6122b7e0a414debcd5ded135ab20'),'filtered_LFP')
catch
	filtered_LFP = LFP;
	for i = 1:width(LFP)
		filtered_LFP(:,i) = 10*bandPass(LFP(:,i),1e4,Inf);
	end
end

% differentiate the LFP
dLFP = LFP;
for i = 1:length(paradigm)
	% first high pass them to remove spikes
	dLFP(:,i) = bandPass(LFP(:,i),Inf,10);
	dLFP(:,i) = -1e4*filtfilt(ones(10,1),10,[0; diff(dLFP(:,i))]);
end

% define limits on data
a = 10e3; z = 50e3;

% extract filters and compute gains
[K1,K1p,K1_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

% we compute the firing machinery gain on the derivative of the LFP to avoid problems with filtering 
[K2,K2p,K2_gain] = extractFilters(dLFP,fA,'use_cache',true,'a',a,'z',z);
[K3,K3p,K3_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);
ft = 1e-3*(1:length(K1)) - .1;

% show transduction i/o curves
ss = 100;
c = parula(max(paradigm)+1);
mean_stim = nanmean(PID(a:z,:));
ms = [min(mean_stim) max(mean_stim)];
ax(2) = subplot(2,5,2); hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(filtered_LFP(a:z,paradigm == i),2);
	x = nanmean(K1p(a:z,paradigm == i),2);
	x = x - nanmean(x);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('Projected Stimulus (V)')
ylabel('ab3 \DeltaLFP (mV)')

% show transduction gain
ax(3) = subplot(2,5,3); hold on
for i = 1:length(paradigm)
	plot(mean_stim(i),K1_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel('\mu_{Odor} (V)')
ylabel('Transduction Gain (mV/V)')
set(gca,'XScale','log','YScale','log')
ff = fit(mean_stim(:),K1_gain(:),'power1','Upper',[Inf -1],'Lower',[0 -1]);
plot(ms,ff(ms),'r')
set(gca,'XScale','log','YScale','log','YLim',[1 100],'XLim',[.1 2])

th = text(.8, 10,'$\sim 1/s$','interpreter','latex','Color','r','FontSize',20);


% show firing machinery I/o curves
ax(7) = subplot(2,5,7); hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = .1*nanmean(K2p(a:z,paradigm == i),2);
	x = x - nanmean(x);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('Projected LFP (mV)')
ylabel('ab3A Firing Rate (Hz)')

% show firing machinery gain
ax(8) = subplot(2,5,8); hold on
for i = 1:length(paradigm)
	plot(mean_stim(i),10*K2_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel('\mu_{Odor} (V)')
ylabel('ab3A Firing Gain (Hz/mV)')
set(gca,'XScale','log','YScale','log','YLim',[1 1e2],'XLim',[.1 2])

clear msg_data
msg_data.PID = PID;
msg_data.LFP = dLFP;
msg_data.fA = fA;
msg_data.paradigm = paradigm;


%  ######   #######  ##    ## ######## ########     ###     ######  ######## 
% ##    ## ##     ## ###   ##    ##    ##     ##   ## ##   ##    ##    ##    
% ##       ##     ## ####  ##    ##    ##     ##  ##   ##  ##          ##    
% ##       ##     ## ## ## ##    ##    ########  ##     ##  ######     ##    
% ##       ##     ## ##  ####    ##    ##   ##   #########       ##    ##    
% ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##    ##    ##    
%  ######   #######  ##    ##    ##    ##     ## ##     ##  ######     ##    


[PID, LFP, fA, ~, orn] = consolidateData(dm.getPath('7955d1ed77512dfe3452b39d71a50e1b'),1);

global_start = 40e3; % 40 seconds
global_end = length(PID) - 5e3; 

% filter the LFP
for i = 1:width(LFP)
	LFP(:,i) = LFP(:,i) - fastFiltFilt(ones(1e4,1),1e4,LFP(:,i));
	LFP(:,i) = LFP(:,i)*10; % to get the units right, now in mV
end

% reshape the LFP signals
block_length = 1e4;
reshaped_LFP = LFP(global_start:end-1e4-1,1:width(PID));
reshaped_LFP = reshape(reshaped_LFP,block_length,width(reshaped_LFP)*length(reshaped_LFP)/block_length);

% also reshape the PID
reshaped_PID = PID(global_start:end-1e4-1,1:width(PID));
reshaped_PID = reshape(reshaped_PID,block_length,width(reshaped_PID)*length(reshaped_PID)/block_length);

% reshape the firing rate signals
reshaped_fA = fA(global_start:end-1e4-1,1:width(PID));
reshaped_fA = reshape(reshaped_fA,block_length,width(reshaped_fA)*length(reshaped_fA)/block_length);


% also reshape the orn ID
reshaped_orn = repmat(orn,length(global_start:length(PID)-1e4-1)/block_length,1);
reshaped_orn = reshaped_orn(:);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];
reshaped_orn(rm_this) = [];

% filter to remove spikes
for i = 1:width(reshaped_LFP)
	reshaped_LFP(:,i) = filtfilt(ones(30,1),30,reshaped_LFP(:,i));
end

% extract filters and find gain
a = 1e3; z = 10e3;
K1 = extractFilters(reshaped_PID,reshaped_LFP,'use_cache',true,'a',a,'z',z);
K2 = extractFilters(reshaped_LFP,reshaped_fA,'use_cache',true,'a',a,'z',z);
K3 = extractFilters(reshaped_PID,reshaped_fA,'use_cache',true,'a',a,'z',z);
ft = 1e-3*(1:length(K1)) - .1;

% remove mean from the LFP for each trial
for i = 1:width(reshaped_LFP)
	reshaped_LFP(:,i) =  reshaped_LFP(:,i) - nanmean(reshaped_LFP(:,i));
end

% average filters for display and to project the stimulus
K1 = nanmean(K1,2);
K2 = nanmean(K2,2); 
K3 = nanmean(K3,2);

% project the stimulus
K1p = NaN*reshaped_fA;
K2p = NaN*reshaped_fA;
K3p = NaN*reshaped_fA;
K2K1p = NaN*reshaped_fA;
for i = 1:width(reshaped_fA)
	K1p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K1,ft);
	K2K1p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),K1p(:,i),K2,ft);
	K2p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_LFP(:,i),K2,ft);
	K3p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K3,ft);
end



% compute the r2 for each trial
r2_K1p = NaN(width(reshaped_LFP),1);
r2_K2p = NaN(width(reshaped_LFP),1);
r2_K2K1p = NaN(width(reshaped_LFP),1);
r2_K3p = NaN(width(reshaped_LFP),1);
for i = 1:length(r2_K1p)
	r2_K1p(i) = rsquare(K1p(1e3:9e3,i),reshaped_LFP(1e3:9e3,i));
	r2_K2p(i) = rsquare(K2p(1e3:9e3,i),reshaped_fA(1e3:9e3,i));
	r2_K2K1p(i) = rsquare(K2K1p(1e3:9e3,i),reshaped_fA(1e3:9e3,i));
	r2_K3p(i) = rsquare(K3p(1e3:9e3,i),reshaped_fA(1e3:9e3,i));
end

ss = 1;
min_r2 = .8; 

% compute gains per trial
lo_gain_LFP = NaN(width(reshaped_PID),1);
hi_gain_LFP = NaN(width(reshaped_PID),1);
for i = 1:width(reshaped_PID)
	y = reshaped_LFP(1e3:4e3,i);
	x = K1p(1e3:4e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		hi_gain_LFP(i) = ff.p1;
	catch
	end

	y = reshaped_LFP(6e3:9e3,i);
	x = K1p(6e3:9e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		lo_gain_LFP(i) = ff.p1;
	catch
	end
end

% compute firing rate gains using Hill functions
if ~exist('lo_gain_firing','var')
	lo_gain_firing = NaN(width(reshaped_PID),1);
	hi_gain_firing = NaN(width(reshaped_PID),1);
	for i = 1:width(reshaped_PID)
		if ~being_published
			textbar(i,width(reshaped_PID))
		end
		y = reshaped_fA(1e3:4e3,i); s = nanmax(y); 	y = y/s;
		x = K2p(1e3:4e3,i);
		try
			ft = fittype('hill2(x,k,n,x_offset)');
			ff = fit(x(:),y(:),ft,'StartPoint',[nanmax(x)/2 2 nanmean(x)],'Lower',[0 1 -Inf],'Upper',[nanmax(x) 10 nanmax(x)],'MaxIter',1e4);
			hi_gain_firing(i) = s*differentiate(ff,ff.k + ff.x_offset);
		catch
		end

		y = reshaped_fA(6e3:9e3,i); s = nanmax(y); 	y = y/s;
		x = K2p(6e3:9e3,i);
		try
			ft = fittype('hill2(x,k,n,x_offset)');
			ff = fit(x(:),y(:),ft,'StartPoint',[nanmax(x)/2 2 nanmean(x)],'Lower',[0 1 -Inf],'Upper',[nanmax(x) 10 nanmax(x)],'MaxIter',1e4);
			lo_gain_firing(i) = s*differentiate(ff,ff.k + ff.x_offset);
		catch
		end
	end
end

% plot transduction i/o curves
ax(4) = subplot(2,5,4); hold on
x = K1p(1e3:4e3,:);
y = reshaped_LFP(1e3:4e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K1p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
for i = 1:width(x)
	x(:,i) = x(:,i) - nanmean(x(:,i));
	y(:,i) = y(:,i) - nanmean(y(:,i));
end
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[1 0 0],'proportional_bins',true,'show_error',false);
x = K1p(6e3:9e3,:);
y = reshaped_LFP(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K1p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
for i = 1:width(x)
	x(:,i) = x(:,i) - nanmean(x(:,i));
	y(:,i) = y(:,i) - nanmean(y(:,i));
end
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[0 0 1],'proportional_bins',true,'show_error',false);
xlabel('Projected Stimulus (V)')
ylabel('ab3 \DeltaLFP (mV)')

% plot transduction gain change
ax(5) = subplot(2,5,5); hold on
x = std(reshaped_PID(1e3:4e3,r2_K1p>.8));
y = hi_gain_LFP(r2_K1p>.8);
plot(ax(5),x,y,'+','Color',[1 opacity opacity])
errorbar(ax(5),nanmean(x),nanmean(y),nanstd(y),'r','LineWidth',4,'Marker','o','MarkerSize',10);
x = std(reshaped_PID(6e3:9e3,r2_K1p>.8));
y = lo_gain_LFP(r2_K1p>.8);
plot(ax(5),x,y,'+','Color',[opacity opacity 1])
errorbar(ax(5),nanmean(x),nanmean(y),nanstd(y),'b','LineWidth',4,'Marker','o','MarkerSize',10);
xlabel('\sigma_{Odor} (V)')
ylabel('Transduction Gain (mV/V)')

% show firing i/o curves
ax(9) = subplot(2,5,9); hold on
x = K2p(1e3:4e3,:);
y = reshaped_fA(1e3:4e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K2p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
for i = 1:width(x)
	x(:,i) = x(:,i) - nanmean(x(:,i));
end
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
[~,data_hi] = plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[1 0 0],'proportional_bins',true,'show_error',false);
x = K2p(6e3:9e3,:);
y = reshaped_fA(6e3:9e3,:);
rm_this = (isnan(sum(y)) | isnan(sum(x)) | (r2_K2p < min_r2)');
x(:,rm_this) = []; y(:,rm_this) = []; 
for i = 1:width(x)
	x(:,i) = x(:,i) - nanmean(x(:,i));
end
x = x(:); y = y(:);
x = x(1:ss:end); y = y(1:ss:end);
[~,data_lo] = plotPieceWiseLinear(x(:),y(:),'nbins',50,'Color',[0 0 1],'proportional_bins',true,'show_error',false);
xlabel('Projected LFP (mV)')
ylabel('ab3A Firing Rate (Hz)')

% plot firing gain change
ax(10) = subplot(2,5,10); hold on; 
x = std(reshaped_PID(1e3:4e3,r2_K2p>.8));
y = hi_gain_firing(r2_K2p>.8);
plot(ax(10),x,y,'+','Color',[1 opacity opacity])
errorbar(ax(10),nanmean(x),nanmean(y),nanstd(y),'r','LineWidth',4,'Marker','o','MarkerSize',10);
x = std(reshaped_PID(6e3:9e3,r2_K2p>.8));
y = lo_gain_firing(r2_K2p>.8);
plot(ax(10),x,y,'+','Color',[opacity opacity 1])
errorbar(ax(10),nanmean(x),nanmean(y),nanstd(y),'b','LineWidth',4,'Marker','o','MarkerSize',10);
xlabel('\sigma_{Odor} (V)')
ylabel('ab3A Firing Gain (Hz/mV)')

% add an explanatory graphic
ax(1) = subplot(1,5,1); hold on
o = imread('../images/fig-LFP-cartoon.png');
imagesc(o);
axis ij
axis image
axis off

prettyFig('fs',16)

% plot example traces of the LFP, firing rate and the stimulus
inset1 = axes;
inset1.Position = [0.13 0.76 0.1 0.1];
plot(inset1,1e-3*(1:5e3),reshaped_PID(1:5e3,1),'k','LineWidth',1.5);

inset2 = axes;
inset2.Position = [0.13 0.54 0.1 0.1];
plot(inset2,1e-3*(1:5e3),reshaped_LFP(1:5e3,1),'k','LineWidth',1.5);

inset3 = axes;
inset3.Position = [0.13 0.24 0.1 0.1];
plot(inset3,1e-3*(1:5e3),reshaped_fA(1:5e3,1),'k','LineWidth',1.5);
set(inset1,'LineWidth',1.5,'box','off','XLim',[0 5])
set(inset2,'LineWidth',1.5,'box','off','XLim',[0 5])
set(inset3,'LineWidth',1.5,'box','off','XLim',[0 5])
xlabel(inset3,'Time (s)')

% make all the plots square
for i = 2:length(ax)
	try
		axis(ax(i),'square')
	catch
	end
end

% move some plots around
ax(1).Position(1) = .01;

ax(5).Position(1) = .85;
ax(4).Position(1) = .68;
ax(9).Position(1) = .68;
ax(10).Position(1) = .85;

ax(2).Position(2) = .53;
ax(3).Position(2) = .53;
ax(4).Position(2) = .53;
ax(5).Position(2) = .53;

th.FontSize = 20;

% fake a axes that covers the entire figure
a = axes;
a.Position = [0 0 1 1];
uistack(a,'bottom')
a.TickLength = [0 0];
a.XLim = [0 1];
a.YLim = [0 1];

% draw a couple of rectangles to box things in
r1 = rectangle('Position',[0.1 0.1 .4 .4],'Curvature',0);
r1.Position = [.255 .025 .34 .95];
r1.FaceColor = [.9 .9 .9];
r1.LineStyle = 'none';
r1.Parent = a;

r2 = rectangle('Position',[0.1 0.1 .4 .4],'Curvature',0);
r2.Position = [.64 .025 .345 .95];
r2.FaceColor = [.9 .9 .9];
r2.LineStyle = 'none';
r2.Parent = a;

% add some text
t1 = text;
t1.String = 'Changing Stimulus Mean';
t1.FontSize = 20;
t1.FontWeight=  'bold';
t1.Parent = a;
t1.Position = [0.3500 0.9100 0];

t2 = text;
t2.String = 'Changing Stimulus Variance';
t2.FontSize = 20;
t2.FontWeight=  'bold';
t2.Parent = a;
t2.Position = [0.7400 0.9100 0];



if being_published	
	snapnow	
	delete(gcf)
end


%% Supplementary Figure
% In this supplementary figure, we show the same result by estimating gain as ratios of standard deviations of output to input

figure('outerposition',[0 0 800 900],'PaperUnits','points','PaperSize',[800 900]); hold on
subplot(2,2,1), hold on
x = mean(msg_data.PID(20e3:45e3,:));
y = std(msg_data.LFP(20e3:45e3,:))./std(msg_data.PID(20e3:45e3,:));
c = parula(length(unique(msg_data.paradigm))+1);
for i = 1:length(x)
	plot(x(i),y(i),'+','Color',c(msg_data.paradigm(i),:));
end
ff = fit(x(:),y(:),'power1','Upper',[Inf -1],'Lower',[0 -1]);
plot(x,ff(x),'r')
set(gca,'XScale','log','YScale','log','YLim',[10 1e3],'XTick',[.1 1],'XLim',[.1 2])
xlabel('\mu_{Odor} (V)')
ylabel('\sigma_{LFP}/\sigma_{Stimulus} (mV/V)')

subplot(2,2,3), hold on
x = mean(msg_data.PID(20e3:45e3,:));
y = 10*std(msg_data.fA(20e3:45e3,:))./std(msg_data.LFP(20e3:45e3,:));
c = parula(length(unique(msg_data.paradigm))+1);
for i = 1:length(x)
	plot(x(i),y(i),'+','Color',c(msg_data.paradigm(i),:));
end
set(gca,'XScale','log','YScale','log','YLim',[1 100],'XTick',[.1 1],'XLim',[.1 2])
ylabel('\sigma_{Firing Rate}/\sigma_{LFP} (Hz/mV)')
xlabel('\mu_{Odor} (V)')

% now do the variance switching experiment
subplot(2,2,2), hold on
y = std(reshaped_LFP(1e3:5e3,:))./std(reshaped_PID(1e3:5e3,:));
x = std(reshaped_PID(1e3:5e3,:));
plot(x,y,'+','Color',[1 opacity opacity])
errorbar(nanmean(x),nanmean(y),nanstd(y),'r','LineWidth',4,'Marker','o','MarkerSize',10);

y = std(reshaped_LFP(6e3:9e3,:))./std(reshaped_PID(6e3:9e3,:));
x = std(reshaped_PID(6e3:9e3,:));
plot(x,y,'+','Color',[opacity opacity 1])
errorbar(nanmean(x),nanmean(y),nanstd(y),'b','LineWidth',4,'Marker','o','MarkerSize',10);
xlabel('\sigma_{Odor} (V)')
ylabel('\sigma_{LFP}/\sigma_{Stimulus} (mV/V)')
set(gca,'YLim',[0 20])

subplot(2,2,4), hold on
y = std(reshaped_fA(1e3:5e3,:))./std(reshaped_LFP(1e3:5e3,:));
x = std(reshaped_PID(1e3:5e3,:));
plot(x,y,'+','Color',[1 opacity opacity])
errorbar(nanmean(x),nanmean(y),nanstd(y),'r','LineWidth',4,'Marker','o','MarkerSize',10);

y = std(reshaped_fA(6e3:9e3,:))./std(reshaped_LFP(6e3:9e3,:));
x = std(reshaped_PID(6e3:9e3,:));
plot(x,y,'+','Color',[opacity opacity 1])
errorbar(nanmean(x),nanmean(y),nanstd(y),'b','LineWidth',4,'Marker','o','MarkerSize',10);
set(gca,'YLim',[0 60])
ylabel('\sigma_{Firing Rate}/\sigma_{LFP} (Hz/mV)')
xlabel('\sigma_{Odor} (V)')

prettyFig
if being_published	
	snapnow	
	delete(gcf)
end

%% Version Info
%
pFooter;


