%% fig_LFP.m
% makes figure showing how gain changes at the LFP level, and from the LFP to the firing machinery
% 


pHeader;

% uses dataManager for data integrity
dm = dataManager;

opacity = .5;

figure('outerposition',[0 0 1600 720],'PaperUnits','points','PaperSize',[1600 720]); hold on

clear cdata
cdata = consolidateData2(dm.getPath('93ba5d68174e3df9f462a1fc48c581da'));
cdata = cleanMSGdata(cdata);

v2struct(cdata)

% [PID, LFP, fA, paradigm,~, ~, AllControlParadigms] = consolidateData(dm.getPath('93ba5d68174e3df9f462a1fc48c581da'),1);
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
[~,K1p,K1_gain] = extractFilters(PID,filtered_LFP,'use_cache',true,'a',a,'z',z);

% we compute the firing machinery gain on the derivative of the LFP to avoid problems with filtering 
[~,K2p,K2_gain] = extractFilters(dLFP,fA,'use_cache',true,'a',a,'z',z);
% [K3,K3p,K3_gain] = extractFilters(PID,fA,'use_cache',true,'a',a,'z',z);
% ft = 1e-3*(1:length(K1)) - .1;

% show transduction i/o curves
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
xlabel('Projected stimulus (V)')
ylabel('ab3 \DeltaLFP (mV)')

% show transduction gain
ax(3) = subplot(2,5,3); hold on
for i = 1:length(paradigm)
	plot(mean_stim(i),K1_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel('\mu_{Stimulus} (V)')
ylabel('ab3 transduction gain (mV/V)')
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
ylabel('ab3A Firing rate (Hz)')

% show firing machinery gain
ax(8) = subplot(2,5,8); hold on
for i = 1:length(paradigm)
	plot(mean_stim(i),10*K2_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel('\mu_{Stimulus} (V)')
ylabel('ab3A Firing gain (Hz/mV)')
set(ax(8),'XScale','log','YScale','log','YLim',[1e-1 1e1],'XLim',[.1 2])

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

clearvars -except dm msg_data being_published opacity ax lo_gain_firing hi_gain_firing lo_gain_firing_corrected hi_gain_firing_corrected lo_gain_total hi_gain_total hi_gain_total_corrected lo_gain_total_corrected

[PID, LFP, fA, ~, orn,~,~,~,~,all_spikes] = consolidateData(dm.getPath('e30707e8e8ef6c0d832eee31eaa585aa'),1);

% compress all_spikes into a 1ms time step
all_spikes = all_spikes';
all_spikes2 = 0*LFP;
for i = 1:width(all_spikes)
	temp = ones(10,1);
	temp2 = filter(temp,10,full(all_spikes(:,i)));
	temp2(temp2>0) = 1;
	temp2(temp2<1) = 0;
	all_spikes2(:,i) = temp2(1:10:end);
end
clear all_spikes

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

% reshape the raw spike times
reshaped_spikes = all_spikes2(global_start:end-1e4-1,1:width(PID));
reshaped_spikes = reshape(reshaped_spikes,block_length,width(reshaped_spikes)*length(reshaped_spikes)/block_length);


% also reshape the orn ID
reshaped_orn = repmat(orn,length(global_start:length(PID)-1e4-1)/block_length,1);
reshaped_orn = reshaped_orn(:);

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];
reshaped_fA(:,rm_this) = [];
reshaped_orn(rm_this) = [];
reshaped_spikes(:,rm_this) = [];

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

for i = 1:width(reshaped_fA)
	K1p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K1,ft);
	K2p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_LFP(:,i),K2,ft);
	K3p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),K3,ft);
end

% correct projections by the mean stimulus in each trial
K1p_corrected = NaN*K1p;
%K2p_corrected = NaN*K2p; % we don't want to correct this -- makes no sense
K3p_corrected = NaN*K3p;
for i = 1:width(reshaped_fA)
	K1p_corrected(1:5e3,i) = K1p(1:5e3,i)/mean(reshaped_PID(1e3:4e3,i));
	K1p_corrected(5e3+1:end,i) = K1p(5e3+1:end,i)/mean(reshaped_PID(6e3:9e3,i));

	%K2p_corrected(1:5e3,i) = K2p(1:5e3,i)/mean(reshaped_PID(1e3:4e3,i));
	%K2p_corrected(5e3+1:end,i) = K2p(5e3+1:end,i)/mean(reshaped_PID(6e3:9e3,i));

	K3p_corrected(1:5e3,i) = K3p(1:5e3,i)/mean(reshaped_PID(1e3:4e3,i));
	K3p_corrected(5e3+1:end,i) = K3p(5e3+1:end,i)/mean(reshaped_PID(6e3:9e3,i));
end

K1p_corrected = K1p_corrected*mean(reshaped_PID(:)); % overall correction to get the units right
%K2p_corrected = K2p_corrected*mean(reshaped_PID(:)); % overall correction to get the units right
K3p_corrected = K3p_corrected*mean(reshaped_PID(:)); % overall correction to get the units right


% compute the r2 for each trial
r2_K1p = NaN(width(reshaped_LFP),1);
r2_K2p = NaN(width(reshaped_LFP),1);
r2_K3p = NaN(width(reshaped_LFP),1);
for i = 1:length(r2_K1p)
	r2_K1p(i) = rsquare(K1p(1e3:9e3,i),reshaped_LFP(1e3:9e3,i));
	r2_K2p(i) = rsquare(K2p(1e3:9e3,i),reshaped_fA(1e3:9e3,i));
	r2_K3p(i) = rsquare(K3p(1e3:9e3,i),reshaped_fA(1e3:9e3,i));
end

ss = 1;
min_r2 = .8; 
ok = (r2_K3p > min_r2 & r2_K1p > min_r2 & r2_K2p > min_r2);

% compute gains per trial on both the corrected and uncorrected LFP predictions
lo_gain_LFP = NaN(width(reshaped_PID),1);
hi_gain_LFP = NaN(width(reshaped_PID),1);
lo_gain_LFP_corrected = NaN(width(reshaped_PID),1);
hi_gain_LFP_corrected = NaN(width(reshaped_PID),1);
for i = 1:width(reshaped_PID)
	y = reshaped_LFP(1e3:4e3,i);
	xc = K1p_corrected(1e3:4e3,i);
	x = K1p(1e3:4e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		hi_gain_LFP(i) = ff.p1;
		ff = fit(xc(:),y(:),'poly1');
		hi_gain_LFP_corrected(i) = ff.p1;
	catch
	end

	y = reshaped_LFP(6e3:9e3,i);
	xc = K1p_corrected(6e3:9e3,i);
	x = K1p(6e3:9e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		lo_gain_LFP(i) = ff.p1;
		ff = fit(xc(:),y(:),'poly1');
		lo_gain_LFP_corrected(i) = ff.p1;
	catch
	end
	clear x y xc
end

% compute firing rate gains using Hill functions
if ~exist('lo_gain_firing','var')
	ft = fittype('hill2(x,k,n,x_offset)');

	lo_gain_firing = NaN(width(reshaped_PID),1);
	hi_gain_firing = NaN(width(reshaped_PID),1);
	lo_gain_total = NaN(width(reshaped_PID),1);
	hi_gain_total = NaN(width(reshaped_PID),1);
	lo_gain_total_corrected = NaN(width(reshaped_PID),1);
	hi_gain_total_corrected = NaN(width(reshaped_PID),1);

	for i = 1:width(reshaped_PID)
		if ~being_published
			textbar(i,width(reshaped_PID))
		end

		% high variance
		y = reshaped_fA(1e3:4e3,i); s = nanmax(y); 	y = y/s;
		x = K2p(1e3:4e3,i);
		%xc = K2p_corrected(1e3:4e3,i);
		x3 = K3p(1e3:4e3,i);
		x3c = K3p_corrected(1e3:4e3,i);

		try
			% LFP ➔ firing 
			ff = fit(x(:),y(:),ft,'StartPoint',[nanmax(x)/2 2 nanmean(x)],'Lower',[0 1 -Inf],'Upper',[nanmax(x) 10 nanmax(x)],'MaxIter',1e4);
			hi_gain_firing(i) = s*differentiate(ff,ff.k + ff.x_offset);
			% ff = fit(xc(:),y(:),ft,'StartPoint',[nanmax(xc)/2 2 nanmean(xc)],'Lower',[0 1 -Inf],'Upper',[nanmax(xc) 10 nanmax(xc)],'MaxIter',1e4);
			% hi_gain_firing_corrected(i) = s*differentiate(ff,ff.k + ff.x_offset);
			% stimulus ➔ firing
			ff = fit(x3(:),y(:),ft,'StartPoint',[nanmax(x)/2 2 nanmean(x)],'Lower',[0 1 -Inf],'Upper',[nanmax(x) 10 nanmax(x)],'MaxIter',1e4);
			hi_gain_total(i) = s*differentiate(ff,ff.k + ff.x_offset);
			ff = fit(x3c(:),y(:),ft,'StartPoint',[nanmax(x3c)/2 2 nanmean(x3c)],'Lower',[0 1 -Inf],'Upper',[nanmax(x3c) 10 nanmax(x3c)],'MaxIter',1e4);
			hi_gain_total_corrected(i) = s*differentiate(ff,ff.k + ff.x_offset);
		catch
		end
		clear x y x3 x3c

		% low variance
		y = reshaped_fA(6e3:9e3,i); s = nanmax(y); 	y = y/s;
		x = K2p(6e3:9e3,i);
		%xc = K2p_corrected(6e3:9e3,i);
		x3 = K3p(6e3:9e3,i);
		x3c = K3p_corrected(6e3:9e3,i);

		try
			% LFP ➔ firing
			ff = fit(x(:),y(:),ft,'StartPoint',[nanmax(x)/2 2 nanmean(x)],'Lower',[0 1 -Inf],'Upper',[nanmax(x) 10 nanmax(x)],'MaxIter',1e4);
			lo_gain_firing(i) = s*differentiate(ff,ff.k + ff.x_offset);
			% ff = fit(xc(:),y(:),ft,'StartPoint',[nanmax(xc)/2 2 nanmean(xc)],'Lower',[0 1 -Inf],'Upper',[nanmax(xc) 10 nanmax(x)],'MaxIter',1e4);
			% lo_gain_firing_corrected(i) = s*differentiate(ff,ff.k + ff.x_offset);
			% stimulus ➔ firing
			ff = fit(x3(:),y(:),ft,'StartPoint',[nanmax(x)/2 2 nanmean(x)],'Lower',[0 1 -Inf],'Upper',[nanmax(x) 10 nanmax(x)],'MaxIter',1e4);
			lo_gain_total(i) = s*differentiate(ff,ff.k + ff.x_offset);
			ff = fit(x3c(:),y(:),ft,'StartPoint',[nanmax(x3c)/2 2 nanmean(x3c)],'Lower',[0 1 -Inf],'Upper',[nanmax(x3c) 10 nanmax(x)],'MaxIter',1e4);
			lo_gain_total_corrected(i) = s*differentiate(ff,ff.k + ff.x_offset);
		catch
		end
		clear x y x3c x3
	end
end

% plot transduction i/o curves
ax(4) = subplot(2,5,4); hold on
x = K1p_corrected(1e3:4e3,:);
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
x = K1p_corrected(6e3:9e3,:);
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
xlabel('Projected stimulus (V)')
ylabel('ab3 \DeltaLFP (mV)')

% plot transduction gain change
ax(5) = subplot(2,5,5); hold on
x = std(reshaped_PID(1e3:4e3,r2_K1p>.8));
y = hi_gain_LFP_corrected(r2_K1p>.8);
plot(ax(5),x,y,'+','Color',[1 opacity opacity])
errorbar(ax(5),nanmean(x),nanmean(y),nanstd(y),'r','LineWidth',4,'Marker','o','MarkerSize',10);
x = std(reshaped_PID(6e3:9e3,r2_K1p>.8));
y = lo_gain_LFP_corrected(r2_K1p>.8);
plot(ax(5),x,y,'+','Color',[opacity opacity 1])
errorbar(ax(5),nanmean(x),nanmean(y),nanstd(y),'b','LineWidth',4,'Marker','o','MarkerSize',10);
xlabel('\sigma_{Stimulus} (V)')
ylabel('ab3 transduction Gain (mV/V)')
set(ax(5),'XLim',[0 .2],'YLim',[0 15])

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
ylabel('ab3A Firing rate (Hz)')

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
xlabel('\sigma_{Stimulus} (V)')
ylabel('ab3A Firing gain (Hz/mV)')
set(ax(10),'XLim',[0 .2],'YLim',[0 50])


prettyFig('fs',16)

% plot example traces of the LFP, firing rate and the stimulus
inset1 = axes;
inset1.Position = [0.13 0.75 0.1 0.1];
plot(inset1,1e-3*(1:5e3),reshaped_PID(1:5e3,1),'k','LineWidth',1.5);

inset2 = axes;
inset2.Position = [0.13 0.45 0.1 0.1];
plot(inset2,1e-3*(1:5e3),reshaped_LFP(1:5e3,1),'k','LineWidth',1.5);

inset3 = axes;
inset3.Position = [0.13 0.12 0.1 0.1];
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
% ax(1).Position(1) = .01;

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

% add some text
t1 = text;
t1.String = 'Changing stimulus mean';
t1.FontSize = 20;
t1.FontWeight=  'bold';
t1.Parent = a;
t1.Position = [0.3500 0.9100 0];

t2 = text;
t2.String = 'Changing stimulus variance';
t2.FontSize = 20;
t2.FontWeight=  'bold';
t2.Parent = a;
t2.Position = [0.7400 0.9100 0];

if being_published	
	snapnow	
	delete(gcf)
end


%% Supplementary Figure
% Supplementary figure for figure 4. 

figure('outerposition',[0 0 800 910],'PaperUnits','points','PaperSize',[800 910]); hold on

% changing mean stimulus, estimating gain directly without a LN model ~~~~~~~~~~~~~~~~~~~~ 
subplot(4,3,1); hold on
title(['Changing stimulus mean' char(10) 'direct gain estimation'])
x = mean(msg_data.PID(20e3:45e3,:));
y = std(msg_data.LFP(20e3:45e3,:))./std(msg_data.PID(20e3:45e3,:));
c = parula(length(unique(msg_data.paradigm))+1);
for i = 1:length(x)
	plot(x(i),.1*y(i),'+','Color',c(msg_data.paradigm(i),:)); % unit correction
end
ff = fit(x(:),.1*y(:),'power1','Upper',[Inf -1],'Lower',[0 -1]);
plot(x,ff(x),'r')
set(gca,'XScale','log','YScale','log','YLim',[10 1e3],'XTick',[.1 1],'XLim',[.1 2])
xlabel('\mu_{Stimulus} (V)')
ylabel('\sigma_{LFP}/\sigma_{Stimulus} (mV/V)')

subplot(4,3,4); hold on
x = mean(msg_data.PID(20e3:45e3,:));
y = 10*std(msg_data.fA(20e3:45e3,:))./std(msg_data.LFP(20e3:45e3,:));
c = parula(length(unique(msg_data.paradigm))+1);
for i = 1:length(x)
	plot(x(i),y(i),'+','Color',c(msg_data.paradigm(i),:));
end
set(gca,'XScale','log','YScale','log','YLim',[.1 10],'XTick',[.1 1],'XLim',[.1 2])
ylabel('\sigma_{Firing rate}/\sigma_{LFP} (Hz/mV)')
xlabel('\mu_{Stimulus} (V)')

%% changing stimulus variance, comparing fold change in gain at LFP to fold change in gain at transduction  ~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~~~~~~~ ~~~~~~

clear ax gains
gains.gL_hi = hi_gain_LFP_corrected(ok);
gains.gL_lo = lo_gain_LFP_corrected(ok);
gains.gF_hi = hi_gain_firing(ok);
gains.gF_lo = lo_gain_firing(ok);
gains.gT_hi = hi_gain_total_corrected(ok);
gains.gT_lo = lo_gain_total_corrected(ok);
gains.s_lo = std(reshaped_PID(6e3:9e3,ok));
gains.s_hi = std(reshaped_PID(1e3:5e3,ok));

ax(5) = subplot(4,3,7); hold on
ax(6) = subplot(4,3,10); hold on
make_plot = false(6,1); make_plot(5:6) = true;
ax = compareGainsInFig4(gains,ax,make_plot);

% move pie chart, prettify it
ax(7).Position = [0.27 0.18 0.07 0.15];
ax(7).Children(1).String = strrep(ax(7).Children(1).String,'%','');
ax(7).Children(1).Position = [.6 -.24];
ax(7).Children(1).Color = 'w';
ax(7).Children(3).String = strrep(ax(7).Children(3).String,'%','');
ax(7).Children(3).Position = [-.55 .21];
ax(7).Children(3).Color = 'w';

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ changing stimulus variance, directly estimating gains and then correcting for change in mean. 

clear gains
gains.gL_hi = std(reshaped_LFP(1e3:5e3,ok))./std(reshaped_PID(1e3:5e3,ok));
gains.gL_lo = std(reshaped_LFP(6e3:9e3,ok))./std(reshaped_PID(6e3:9e3,ok));
gains.gF_hi = std(reshaped_fA(1e3:5e3,ok))./std(reshaped_LFP(1e3:5e3,ok));
gains.gF_lo = std(reshaped_fA(6e3:9e3,ok))./std(reshaped_LFP(6e3:9e3,ok));
gains.gT_hi = std(reshaped_fA(1e3:5e3,ok))./std(reshaped_PID(1e3:5e3,ok));
gains.gT_lo = std(reshaped_fA(6e3:9e3,ok))./std(reshaped_PID(6e3:9e3,ok));
gains.s_lo = std(reshaped_PID(6e3:9e3,ok));
gains.s_hi = std(reshaped_PID(1e3:5e3,ok));
% correct for change in mean
gains.gL_hi = gains.gL_hi.*mean(reshaped_PID(1e3:5e3,ok));
gains.gL_lo = gains.gL_lo.*mean(reshaped_PID(6e3:9e3,ok));
gains.gT_hi = gains.gT_hi.*mean(reshaped_PID(1e3:5e3,ok));
gains.gT_lo = gains.gT_lo.*mean(reshaped_PID(6e3:9e3,ok));


subplot(4,3,8); hold on
plot(gains.s_hi,gains.gL_hi,'+','Color',[1 opacity opacity])
errorbar(nanmean(gains.s_hi),nanmean(gains.gL_hi),nanstd(gains.gL_hi),'r','LineWidth',4,'Marker','o','MarkerSize',10);
plot(gains.s_lo,gains.gL_lo,'+','Color',[opacity opacity 1]);
errorbar(nanmean(gains.s_lo),nanmean(gains.gL_lo),nanstd(gains.gL_lo),'b','LineWidth',4,'Marker','o','MarkerSize',10);
xlabel('\sigma_{Stimulus} (V)')
ylabel(['(\sigma_{LFP}/\sigma_{Stimulus})' char(10) '\times \mu_{Stimulus} (a.u.)'])
set(gca,'XLim',[0 0.2],'YLim',[0 5])

subplot(4,3,11); hold on
plot(gains.s_hi,gains.gF_hi,'+','Color',[1 opacity opacity])
errorbar(nanmean(gains.s_hi),nanmean(gains.gF_hi),nanstd(gains.gF_hi),'r','LineWidth',4,'Marker','o','MarkerSize',10);
plot(gains.s_lo,gains.gF_lo,'+','Color',[opacity opacity 1]);
errorbar(nanmean(gains.s_lo),nanmean(gains.gF_lo),nanstd(gains.gF_lo),'b','LineWidth',4,'Marker','o','MarkerSize',10);
xlabel('\sigma_{Stimulus} (V)')
ylabel('\sigma_{Firing}/\sigma_{LFP} (a.u.)')
set(gca,'XLim',[0 0.2],'YLim',[0 50])

clear ax 
ax(5) = subplot(4,3,9); hold on
ax(6) = subplot(4,3,12); hold on
make_plot = false(6,1); make_plot(5:6) = true;
ax = compareGainsInFig4(gains,ax,make_plot);

% move pie chart, prettify it
ax(7).Position = [0.86 0.18 0.07 0.15];
ax(7).Children(1).String = strrep(ax(7).Children(1).String,'%','');
ax(7).Children(1).Position = [.6 -.24];
ax(7).Children(1).Color = 'w';
ax(7).Children(3).String = strrep(ax(7).Children(3).String,'%','');
ax(7).Children(3).Position = [-.55 .21];
ax(7).Children(3).Color = 'w';

%% Weber's Law in transduction for other odors/receptors \
% define what we want to work on
data_hashes = {'bcd4cf4fe12817d084a2b06f981161ee','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1','f11c4a5792d0c9fec7c40fd6aa2fce40'};
odour_names = {'1-pentanol','1-pentanol','2-butanone','isoamyl-acetate'};
orn_names = {'ab3A','ab2A','ab2A','pb1A'};

clear plot_here
plot_here(1) = subplot(4,3,2); hold on
plot_here(2) = subplot(4,3,3); hold on
plot_here(3) = subplot(4,3,5); hold on
plot_here(4) = subplot(4,3,6); hold on

% core loop
for i = length(data_hashes):-1:1
	clear cdata
	cdata = consolidateData2(dm.getPath(data_hashes{i}));
	cdata = cleanMSGdata(cdata);

	% plot gain as we normally calculate it
	clear ph
	ph = plot_here(i);
	plotMSGGain(cdata,ph);

	t = [orn_names{i} char(10) odour_names{i}];
	title(ph,t);
end

prettyFig('fs',14,'lw',1.5);

if being_published
	snapnow
	delete(gcf)
end

% ######## ##     ## ######## ########     ###    
% ##        ##   ##     ##    ##     ##   ## ##   
% ##         ## ##      ##    ##     ##  ##   ##  
% ######      ###       ##    ########  ##     ## 
% ##         ## ##      ##    ##   ##   ######### 
% ##        ##   ##     ##    ##    ##  ##     ## 
% ######## ##     ##    ##    ##     ## ##     ## 

% ########  ##        #######  ########  ######  
% ##     ## ##       ##     ##    ##    ##    ## 
% ##     ## ##       ##     ##    ##    ##       
% ########  ##       ##     ##    ##     ######  
% ##        ##       ##     ##    ##          ## 
% ##        ##       ##     ##    ##    ##    ## 
% ##        ########  #######     ##     ######  


%% Extra plots showing modular gain
% First, we compute gains directly, and don't correct for the change in mean

clear gains
gains.gL_hi = std(reshaped_LFP(1e3:5e3,ok))./std(reshaped_PID(1e3:5e3,ok));
gains.gL_lo = std(reshaped_LFP(6e3:9e3,ok))./std(reshaped_PID(6e3:9e3,ok));
gains.gF_hi = std(reshaped_fA(1e3:5e3,ok))./std(reshaped_LFP(1e3:5e3,ok));
gains.gF_lo = std(reshaped_fA(6e3:9e3,ok))./std(reshaped_LFP(6e3:9e3,ok));
gains.gT_hi = std(reshaped_fA(1e3:5e3,ok))./std(reshaped_PID(1e3:5e3,ok));
gains.gT_lo = std(reshaped_fA(6e3:9e3,ok))./std(reshaped_PID(6e3:9e3,ok));
gains.s_lo = std(reshaped_PID(6e3:9e3,ok));
gains.s_hi = std(reshaped_PID(1e3:5e3,ok));

ax = compareGainsInFig4(gains);
suptitle('Gains directly estimated; no correction for mean')

%%
% Next, we compute gains directly, but correct the LFP and overall gains for the change in the mean

% correct for change in mean
gains.gL_hi = gains.gL_hi.*mean(reshaped_PID(1e3:5e3,ok));
gains.gL_lo = gains.gL_lo.*mean(reshaped_PID(6e3:9e3,ok));
gains.gT_hi = gains.gT_hi.*mean(reshaped_PID(1e3:5e3,ok));
gains.gT_lo = gains.gT_lo.*mean(reshaped_PID(6e3:9e3,ok));

ax = compareGainsInFig4(gains);
suptitle('Gains directly estimated; corrected for change in mean')

%%
% Now, we compute gains using a LN model, and don't correct for the change in mean. 

clear gains
gains.gL_hi = hi_gain_LFP(ok);
gains.gL_lo = lo_gain_LFP(ok);
gains.gF_hi = hi_gain_firing(ok);
gains.gF_lo = lo_gain_firing(ok);
gains.gT_hi = hi_gain_total(ok);
gains.gT_lo = lo_gain_total(ok);
gains.s_lo = std(reshaped_PID(6e3:9e3,ok));
gains.s_hi = std(reshaped_PID(1e3:5e3,ok));

ax = compareGainsInFig4(gains);
suptitle('Gains estimated from LN model; no correction for mean')

%%
% Finally, we compute gains using a LN mode, and correct for the change in the mean.

clear gains
gains.gL_hi = hi_gain_LFP_corrected(ok);
gains.gL_lo = lo_gain_LFP_corrected(ok);
gains.gF_hi = hi_gain_firing(ok);
gains.gF_lo = lo_gain_firing(ok);
gains.gT_hi = hi_gain_total_corrected(ok);
gains.gT_lo = lo_gain_total_corrected(ok);
gains.s_lo = std(reshaped_PID(6e3:9e3,ok));
gains.s_hi = std(reshaped_PID(1e3:5e3,ok));

ax = compareGainsInFig4(gains);
suptitle('Gains estimated from LN model; corrected for change in mean')


%% Version Info
%
pFooter;


