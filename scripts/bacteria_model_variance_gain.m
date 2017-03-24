


pHeader;

%% Does the bacteria model do variance gain control? 
%



% get the variance switching data
clear VSdata
[VSdata.PID, VSdata.LFP, VSdata.fA, VSdata.paradigm, VSdata.orn, VSdata.fly] = consolidateData(getPath(dataManager,'e30707e8e8ef6c0d832eee31eaa585aa'),1);
% remove baseline from stimulus
VSdata.PID = bsxfun(@minus, VSdata.PID, min(VSdata.PID));
VSdata.LFP = bsxfun(@minus, VSdata.LFP, mean(VSdata.LFP(1:4000,:)));

% get the variance switching data filters 
load(getPath(dataManager,'457ee16a326f47992e35a7d5281f9cc4'));
VSdata.K1 = nanmean(K1,2);


PID = VSdata.PID;
LFP = VSdata.LFP;

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

% throw our NaNs globally. so we're throwing out epochs where the data is incomplete
rm_this = isnan(sum(reshaped_LFP));
reshaped_LFP(:,rm_this) = [];
reshaped_PID(:,rm_this) = [];

% filter to remove spikes
for i = 1:width(reshaped_LFP)
	reshaped_LFP(:,i) = filtfilt(ones(30,1),30,reshaped_LFP(:,i));
end

% project using filter
ft = 1e-3*(1:length(VSdata.K1)) - .1;
for i = 1:width(reshaped_PID)
	K1p(:,i) = convolve(1e-3*(1:length(reshaped_PID)),reshaped_PID(:,i),VSdata.K1,ft);
end

% fit a bacteria model to this
clear fd
for i = 1:size(reshaped_PID,2)
	fd(i).stimulus = reshaped_PID(:,i);
	fd(i).response = reshaped_LFP(:,i);
	fd(i).response(1:1e3) = NaN;
end

clear p
p.B = 9.117;
p.e_L = 0.9961;
p.w0 = 18458;
p.K_1 = 0.09933;
p.K_2 = 1.25;
p.K_tau = 82.8;
p.n = 1;
p.output_scale = -33.34;
p.output_offset = 16.503;

% generate response using this model 
model_prediction = bacteriaModelX(reshaped_PID,p);

% compute the i.o curves

x = model_prediction(1e3:5e3,:);
y = reshaped_LFP(1e3:5e3,:);
[~,data_hi] = plotPieceWiseLinear(x,y,'nbins',50,'make_plot',false);

% low contrast
x = model_prediction(6e3:9e3,:);
y = reshaped_LFP(6e3:9e3,:);
[~,data_lo] = plotPieceWiseLinear(x,y,'nbins',50,'make_plot',false);

data_lo.x = data_lo.x - mean(data_lo.x);
data_hi.x = data_hi.x - mean(data_hi.x);
data_lo.y = data_lo.y - mean(data_lo.y);
data_hi.y = data_hi.y - mean(data_hi.y);


figure('outerposition',[0 0 500 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
plot(data_lo.x,data_lo.y,'b')
plot(data_hi.x,data_hi.y,'r')

prettyFig();

if being_published
	snapnow
	delete(gcf)
end

%% Version Info
%
pFooter;

