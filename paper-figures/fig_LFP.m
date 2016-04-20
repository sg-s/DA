


pHeader;


figure('outerposition',[0 0 1700 700],'PaperUnits','points','PaperSize',[1700 700]); hold on


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
subplot(2,5,1), hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(filtered_LFP(a:z,paradigm == i),2);
	x = nanmean(K1p(a:z,paradigm == i),2);
	x = x - nanmean(x);
	% s = nanmean(PID(a:z,paradigm==i),2);
	% x = x + nanmean(s);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('Projected Stimulus (V)')
ylabel('\DeltaLPF (mV)')

% show transduction gain
subplot(2,5,2), hold on
for i = 1:length(paradigm)
	plot(mean_stim(i),K1_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel('\mu_{Stimulus} (V)')
ylabel('Transduction Gain (mV/V)')
set(gca,'XScale','log','YScale','log')
ff = fit(mean_stim(:),K1_gain(:),'power1','Upper',[Inf -1],'Lower',[0 -1]);
plot(ms,ff(ms),'r')
set(gca,'XScale','log','YScale','log','YLim',[1 100],'XLim',[.1 2])


% show firing machinery I/o curves
subplot(2,5,6), hold on
for i = 1:max(paradigm) % iterate over all paradigms 
	y = nanmean(fA(a:z,paradigm == i),2);
	x = .1*nanmean(K2p(a:z,paradigm == i),2);
	x = x - nanmean(x);
	plotPieceWiseLinear(x,y,'nbins',50,'Color',c(i,:));
end
xlabel('Projected LFP (mV)')
ylabel('Firing Rate (Hz)')


subplot(2,5,7), hold on
for i = 1:length(paradigm)
	plot(mean_stim(i),10*K2_gain(i),'+','Color',c(paradigm(i),:))
end
xlabel('\mu_{Stimulus} (V)')
ylabel('Firing Gain (Hz/mV)')
set(gca,'XScale','log','YScale','log','YLim',[1 1e2],'XLim',[.1 2])

%  ######   #######  ##    ## ######## ########     ###     ######  ######## 
% ##    ## ##     ## ###   ##    ##    ##     ##   ## ##   ##    ##    ##    
% ##       ##     ## ####  ##    ##    ##     ##  ##   ##  ##          ##    
% ##       ##     ## ## ## ##    ##    ########  ##     ##  ######     ##    
% ##       ##     ## ##  ####    ##    ##   ##   #########       ##    ##    
% ##    ## ##     ## ##   ###    ##    ##    ##  ##     ## ##    ##    ##    
%  ######   #######  ##    ##    ##    ##     ## ##     ##  ######     ##    


path_name = '/local-data/DA-paper/switching/variance/v2/';
[PID, LFP, fA, ~, orn] = consolidateData(path_name,1);

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
lo_gain_firing = NaN(width(reshaped_PID),1);
hi_gain_firing = NaN(width(reshaped_PID),1);
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

	y = reshaped_fA(1e3:4e3,i);
	x = K2p(1e3:4e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		hi_gain_firing(i) = ff.p1;
	catch
	end

	y = reshaped_fA(6e3:9e3,i);
	x = K2p(6e3:9e3,i);
	try
		ff = fit(x(:),y(:),'poly1');
		lo_gain_firing(i) = ff.p1;
	catch
	end
end

% plot transduction i/o curves
subplot(2,5,4), hold on
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
ylabel('\DeltaLPF (mV)')

% plot transduction gain change
subplot(2,5,5), hold on
x = std(reshaped_PID(1e3:4e3,r2_K1p>.8));
plot(x,hi_gain_LFP(r2_K1p>.8),'r+')
x = std(reshaped_PID(6e3:9e3,r2_K1p>.8));
plot(x,lo_gain_LFP(r2_K1p>.8),'b+')
xlabel('\sigma_{Stimulus} (V)')
ylabel('Firing Gain (Hz/mV)')


% show firing i/o curves
subplot(2,5,9), hold on
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
ylabel('Firing Rate (Hz)')

% plot firing gain change
subplot(2,5,10), hold on
x = std(reshaped_PID(1e3:4e3,r2_K2p>.8));
plot(x,hi_gain_firing(r2_K2p>.8),'r+')
x = std(reshaped_PID(6e3:9e3,r2_K2p>.8));
plot(x,lo_gain_firing(r2_K2p>.8),'b+')
xlabel('\sigma_{Stimulus} (V)')
ylabel('Firing Gain (Hz/mV)')


prettyFig('fs',16)

if being_published	
	snapnow	
	delete(gcf)
end


%% Version Info
%
pFooter;


