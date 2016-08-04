%% fig_LFP.m
% makes figure showing how gain changes at the LFP level, and from the LFP to the firing machinery
% 


pHeader;


opacity = .5;

figure('outerposition',[0 0 1201 900],'PaperUnits','points','PaperSize',[1201 900]); hold on
msg_cartoon_plot = subplot(3,4,1:2); hold on
var_cartoon_plot = subplot(3,4,3:4); hold on

msg_lfp_io_plot = subplot(3,4,5); hold on
msg_lfp_gain_plot = subplot(3,4,6); hold on
msg_firing_io_plot = subplot(3,4,9); hold on
msg_firing_gain_plot = subplot(3,4,10); hold on

var_lfp_io_plot = subplot(3,4,7); hold on
var_lfp_gain_plot = subplot(3,4,8); hold on
var_firing_io_plot = subplot(3,4,11); hold on
var_firing_gain_plot = subplot(3,4,12); hold on

% show the cartoons
axes(msg_cartoon_plot)
o = imread('../images/mean_modular.png');
imagesc(o);
axis ij
axis image
axis off

axes(var_cartoon_plot)
o = imread('../images/variance_modular.png');
imagesc(o);
axis ij
axis image
axis off


clear cdata
cdata = consolidateData2(getPath(dataManager,'93ba5d68174e3df9f462a1fc48c581da'));
cdata = cleanMSGdata(cdata);

v2struct(cdata)

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
	load(getPath(dataManager,'213e6122b7e0a414debcd5ded135ab20'),'filtered_LFP')
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
axes(msg_lfp_io_plot)
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(filtered_LFP(a:z,paradigm == i),2);
	x = nanmean(K1p(a:z,paradigm == i),2);
	x = x - nanmean(x);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel(msg_lfp_io_plot,'Projected stimulus (V)')
ylabel(msg_lfp_io_plot,'ab3 \DeltaLFP (mV)')

% show transduction gain
for i = 1:length(paradigm)
	plot(msg_lfp_gain_plot,mean_stim(i),K1_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel(msg_lfp_gain_plot,'\mu_{Stimulus} (V)')
ylabel(msg_lfp_gain_plot,'ab3 transduction gain (mV/V)')
set(msg_lfp_gain_plot,'XScale','log','YScale','log')
ff = fit(mean_stim(:),K1_gain(:),'power1','Upper',[Inf -1],'Lower',[0 -1]);
plot(msg_lfp_gain_plot,ms,ff(ms),'r')
set(msg_lfp_gain_plot,'XScale','log','YScale','log','YLim',[1 100],'XLim',[.1 2])
axes(msg_lfp_gain_plot)
th = text(.8, 10,'$\sim 1/s$','interpreter','latex','Color','r','FontSize',20);


% show firing machinery I/o curves
axes(msg_firing_io_plot)
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = .1*nanmean(K2p(a:z,paradigm == i),2);
	x = x - nanmean(x);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel(msg_firing_io_plot,'Projected LFP (mV)')
ylabel(msg_firing_io_plot,'ab3A Firing rate (Hz)')

% show firing machinery gain
for i = 1:length(paradigm)
	plot(msg_firing_gain_plot,mean_stim(i),10*K2_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel(msg_firing_gain_plot,'\mu_{Stimulus} (V)')
ylabel(msg_firing_gain_plot,'ab3A Firing gain (Hz/mV)')
set(msg_firing_gain_plot,'XScale','log','YScale','log','YLim',[1e-1 1e1],'XLim',[.1 2])

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

clearvars -except dm msg_data being_published opacity ax lo_gain_firing hi_gain_firing lo_gain_firing_corrected hi_gain_firing_corrected lo_gain_total hi_gain_total hi_gain_total_corrected lo_gain_total_corrected msg_cartoon_plot var_cartoon_plot msg_lfp_io_plot msg_lfp_gain_plot msg_firing_gain_plot var_lfp_io_plot var_lfp_gain_plot var_firing_io_plot var_firing_gain_plot msg_firing_io_plot


[PID, LFP, fA, ~, orn,~,~,~,~,all_spikes] = consolidateData(getPath(dataManager,'e30707e8e8ef6c0d832eee31eaa585aa'),1);

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
axes(var_lfp_io_plot)
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
xlabel(var_lfp_io_plot,'Projected stimulus (V)')
ylabel(var_lfp_io_plot,'ab3 \DeltaLFP (mV)')

% plot transduction gain change
x = std(reshaped_PID(1e3:4e3,r2_K1p>.8));
y = hi_gain_LFP_corrected(r2_K1p>.8);
plot(var_lfp_gain_plot,x,y,'+','Color',[1 opacity opacity])
errorbar(var_lfp_gain_plot,nanmean(x),nanmean(y),nanstd(y),'r','LineWidth',4,'Marker','o','MarkerSize',10);
x = std(reshaped_PID(6e3:9e3,r2_K1p>.8));
y = lo_gain_LFP_corrected(r2_K1p>.8);
plot(var_lfp_gain_plot,x,y,'+','Color',[opacity opacity 1])
errorbar(var_lfp_gain_plot,nanmean(x),nanmean(y),nanstd(y),'b','LineWidth',4,'Marker','o','MarkerSize',10);
xlabel(var_lfp_gain_plot,'\sigma_{Stimulus} (V)')
ylabel(var_lfp_gain_plot,'ab3 transduction Gain (mV/V)')
set(var_lfp_gain_plot,'XLim',[0 .2],'YLim',[0 15])

% show firing i/o curves
axes(var_firing_io_plot)
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
xlabel(var_firing_io_plot,'Projected LFP (mV)')
ylabel(var_firing_io_plot,'ab3A Firing rate (Hz)')

% plot firing gain change
x = std(reshaped_PID(1e3:4e3,r2_K2p>.8));
y = hi_gain_firing(r2_K2p>.8);
plot(var_firing_gain_plot,x,y,'+','Color',[1 opacity opacity])
errorbar(var_firing_gain_plot,nanmean(x),nanmean(y),nanstd(y),'r','LineWidth',4,'Marker','o','MarkerSize',10);
x = std(reshaped_PID(6e3:9e3,r2_K2p>.8));
y = lo_gain_firing(r2_K2p>.8);
plot(var_firing_gain_plot,x,y,'+','Color',[opacity opacity 1])
errorbar(var_firing_gain_plot,nanmean(x),nanmean(y),nanstd(y),'b','LineWidth',4,'Marker','o','MarkerSize',10);
xlabel(var_firing_gain_plot,'\sigma_{Stimulus} (V)')
ylabel(var_firing_gain_plot,'ab3A Firing gain (Hz/mV)')
set(var_firing_gain_plot,'XLim',[0 .2],'YLim',[0 50])


prettyFig('fs',16)
 
% some cosmetic fixes
var_firing_io_plot.Box = 'off';
msg_lfp_io_plot.Position(1) = .1;
msg_firing_io_plot.Position(1) = .1;
msg_firing_gain_plot.Position(1) = .31;
msg_lfp_gain_plot.Position(1) = .31;

var_lfp_io_plot.Position(1) = .57;
var_firing_io_plot.Position(1) = .57;
var_firing_gain_plot.Position(1) = .79;
var_lfp_gain_plot.Position(1) = .79;

msg_cartoon_plot.Position= [.03 .62 .5 .25];
var_cartoon_plot.Position= [.53 .62 .5 .25];

uistack(msg_cartoon_plot,'bottom');
uistack(var_cartoon_plot,'bottom');

if being_published	
	snapnow	
	delete(gcf)
end

%% Supplementary Figure
% Supplementary figure for figure 4. 

figure('outerposition',[0 0 1200 801],'PaperUnits','points','PaperSize',[1200 801]); hold on

% changing mean stimulus, estimating gain directly without a LN model ~~~~~~~~~~~~~~~~~~~~ 
% subplot(4,3,1); hold on
% title(['Changing stimulus mean' char(10) 'direct gain estimation'])
% x = mean(msg_data.PID(20e3:45e3,:));
% y = std(msg_data.LFP(20e3:45e3,:))./std(msg_data.PID(20e3:45e3,:));
% c = parula(length(unique(msg_data.paradigm))+1);
% for i = 1:length(x)
% 	plot(x(i),.1*y(i),'+','Color',c(msg_data.paradigm(i),:)); % unit correction
% end
% ff = fit(x(:),.1*y(:),'power1','Upper',[Inf -1],'Lower',[0 -1]);
% plot(x,ff(x),'r')
% set(gca,'XScale','log','YScale','log','YLim',[10 1e3],'XTick',[.1 1],'XLim',[.1 2])
% xlabel('\mu_{Stimulus} (V)')
% ylabel('\sigma_{LFP}/\sigma_{Stimulus} (mV/V)')

% subplot(4,3,4); hold on
% x = mean(msg_data.PID(20e3:45e3,:));
% y = 10*std(msg_data.fA(20e3:45e3,:))./std(msg_data.LFP(20e3:45e3,:));
% c = parula(length(unique(msg_data.paradigm))+1);
% for i = 1:length(x)
% 	plot(x(i),y(i),'+','Color',c(msg_data.paradigm(i),:));
% end
% set(gca,'XScale','log','YScale','log','YLim',[.1 10],'XTick',[.1 1],'XLim',[.1 2])
% ylabel('\sigma_{Firing rate}/\sigma_{LFP} (Hz/mV)')
% xlabel('\mu_{Stimulus} (V)')

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

ax(5) = subplot(2,3,3); hold on
ax(6) = subplot(2,3,6); hold on
make_plot = false(6,1); make_plot(5:6) = true;
ax = compareGainsInFig4(gains,ax,make_plot);
set(ax(5),'XLim',[0 3],'YLim',[0 3])

% move pie chart, prettify it
%ax(7).Position = [0.27 0.18 0.07 0.15];
ax(7).Children(1).String = strrep(ax(7).Children(1).String,'%','');
ax(7).Children(1).Position = [.6 -.24];
ax(7).Children(1).Color = 'w';
ax(7).Children(3).String = strrep(ax(7).Children(3).String,'%','');
ax(7).Children(3).Position = [-.55 .21];
ax(7).Children(3).Color = 'w';

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ changing stimulus variance, directly estimating gains and then correcting for change in mean. 

% clear gains
% gains.gL_hi = std(reshaped_LFP(1e3:5e3,ok))./std(reshaped_PID(1e3:5e3,ok));
% gains.gL_lo = std(reshaped_LFP(6e3:9e3,ok))./std(reshaped_PID(6e3:9e3,ok));
% gains.gF_hi = std(reshaped_fA(1e3:5e3,ok))./std(reshaped_LFP(1e3:5e3,ok));
% gains.gF_lo = std(reshaped_fA(6e3:9e3,ok))./std(reshaped_LFP(6e3:9e3,ok));
% gains.gT_hi = std(reshaped_fA(1e3:5e3,ok))./std(reshaped_PID(1e3:5e3,ok));
% gains.gT_lo = std(reshaped_fA(6e3:9e3,ok))./std(reshaped_PID(6e3:9e3,ok));
% gains.s_lo = std(reshaped_PID(6e3:9e3,ok));
% gains.s_hi = std(reshaped_PID(1e3:5e3,ok));
% % correct for change in mean
% gains.gL_hi = gains.gL_hi.*mean(reshaped_PID(1e3:5e3,ok));
% gains.gL_lo = gains.gL_lo.*mean(reshaped_PID(6e3:9e3,ok));
% gains.gT_hi = gains.gT_hi.*mean(reshaped_PID(1e3:5e3,ok));
% gains.gT_lo = gains.gT_lo.*mean(reshaped_PID(6e3:9e3,ok));


% subplot(4,3,8); hold on
% plot(gains.s_hi,gains.gL_hi,'+','Color',[1 opacity opacity])
% errorbar(nanmean(gains.s_hi),nanmean(gains.gL_hi),nanstd(gains.gL_hi),'r','LineWidth',4,'Marker','o','MarkerSize',10);
% plot(gains.s_lo,gains.gL_lo,'+','Color',[opacity opacity 1]);
% errorbar(nanmean(gains.s_lo),nanmean(gains.gL_lo),nanstd(gains.gL_lo),'b','LineWidth',4,'Marker','o','MarkerSize',10);
% xlabel('\sigma_{Stimulus} (V)')
% ylabel(['(\sigma_{LFP}/\sigma_{Stimulus})' char(10) '\times \mu_{Stimulus} (a.u.)'])
% set(gca,'XLim',[0 0.2],'YLim',[0 5])

% subplot(4,3,11); hold on
% plot(gains.s_hi,gains.gF_hi,'+','Color',[1 opacity opacity])
% errorbar(nanmean(gains.s_hi),nanmean(gains.gF_hi),nanstd(gains.gF_hi),'r','LineWidth',4,'Marker','o','MarkerSize',10);
% plot(gains.s_lo,gains.gF_lo,'+','Color',[opacity opacity 1]);
% errorbar(nanmean(gains.s_lo),nanmean(gains.gF_lo),nanstd(gains.gF_lo),'b','LineWidth',4,'Marker','o','MarkerSize',10);
% xlabel('\sigma_{Stimulus} (V)')
% ylabel('\sigma_{Firing}/\sigma_{LFP} (a.u.)')
% set(gca,'XLim',[0 0.2],'YLim',[0 50])

% clear ax 
% ax(5) = subplot(4,3,9); hold on
% ax(6) = subplot(4,3,12); hold on
% make_plot = false(6,1); make_plot(5:6) = true;
% ax = compareGainsInFig4(gains,ax,make_plot);

% % move pie chart, prettify it
% ax(7).Position = [0.86 0.18 0.07 0.15];
% ax(7).Children(1).String = strrep(ax(7).Children(1).String,'%','');
% ax(7).Children(1).Position = [.6 -.24];
% ax(7).Children(1).Color = 'w';
% ax(7).Children(3).String = strrep(ax(7).Children(3).String,'%','');
% ax(7).Children(3).Position = [-.55 .21];
% ax(7).Children(3).Color = 'w';

%% Weber's Law in transduction for other odors/receptors \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% define what we want to work on~~~~
data_hashes = {'bcd4cf4fe12817d084a2b06f981161ee','cd6753c0e4cf02895cd5e2c5cb58aa1a','3ea08ccfa892c6545d74bbdaaa6cbee1','f11c4a5792d0c9fec7c40fd6aa2fce40'};
odour_names = {'1-pentanol','1-pentanol','2-butanone','isoamyl-acetate'};
orn_names = {'ab3A','ab2A','ab2A','pb1A'};

clear plot_here
plot_here(1) = subplot(2,3,1); hold on
plot_here(2) = subplot(2,3,2); hold on
plot_here(3) = subplot(2,3,4); hold on
plot_here(4) = subplot(2,3,5); hold on

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

plot_here(4).XLim = [3e-3 3];

prettyFig('fs',14,'lw',1.5,'FixLogX',true);

for i = 1:length(plot_here)
	plot_here(i).XLim(2) = plot_here(i).XLim(2)*1.06;
	deintersectAxes(plot_here(i));
end

deintersectAxes(ax(6))
uistack(ax(7),'top')

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


