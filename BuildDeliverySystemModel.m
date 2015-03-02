% BuildDeliverySystemModel.m
% meant to be run by gaussian constructor.m
% builds a non-parameteric model for the delivery system
% 
% created by Srinivas Gorur-Shandilya at 4:19 , 02 March 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [K_MFC,K_PID,p_hill] = BuildDeliverySystemModel(data,ControlParadigm,use_this)
% core parameters
MFC_Scale = 100;
Total_Flow = 2000;


% subsample data
time = 1e-4*(1:length(data(use_this).PID));
t = time(1:10:end);
PID = zeros(width(data(use_this).PID),length(data(use_this).PID)/10);
MFC = zeros(width(data(use_this).PID),length(data(use_this).PID)/10);
MFC_Control = (ControlParadigm(use_this).Outputs(use_this,1:10:end));
for i = 1:width(data(use_this).PID)
	PID(i,:) = interp1(time,data(use_this).PID(i,:),t);
	MFC(i,:) = interp1(time,data(use_this).MFC500(i,:),t);
end
clear time
time = t;


% extract all MFC filters
K_MFC = zeros(width(PID),600);
for j = 1:width(PID)
	temp = FindBestFilter(MFC_Control,MFC(j,:),[],'regmax=1;','regmin=1;','filter_length=799;','offset=100;');
	K_MFC(j,:) = temp(101:end-100);
end


% make linear predictions
MFC_pred = NaN*MFC;
s =[];
for i = 1:width(MFC)
	MFC_pred(i,:) = filter(K_MFC(i,:),1,MFC_Control);

	% fix trivial scaling
	x = MFC_pred(i,:); x = x(:);
	y = MFC(i,:); y = y(:);
	f = fit(x,y,'poly1');
	MFC_pred(i,:) = f(MFC_pred(i,:));
	s = [s f.p1];
end
s = mean(s);
K_MFC = K_MFC*s;

dil = MFC*100; % 1V = 100mL/min
dil = dil./(dil + 2000);

% extract all PID filters
K_PID = zeros(width(PID),600);
for j = 1:width(PID)
	temp = FindBestFilter(dil(j,:),PID(j,:),[],'regmax=1;','regmin=1;','filter_length=799;','offset=100;');
	K_PID(j,:) = temp(101:end-100);
end


% make linear predictions
PID_pred = NaN*dil;
s =[];
for i = 1:width(MFC)
	PID_pred(i,:) = filter(K_PID(i,:),1,dil(i,:));

	% fix trivial scaling
	x = PID_pred(i,:); x = x(:);
	y = PID(i,:); y = y(:);
	f = fit(x,y,'poly1');
	PID_pred(i,:) = f(PID_pred(i,:));
	s = [s f.p1];
end
s = mean(s);
K_PID = K_PID*s;



% now extract one non-linearity
K_MFC = mean2(K_MFC3(2:end,:));
K_PID = mean2(K_PID(2:end,:));

PID_pred = DeliverySystemModel(MFC_Control);

d.stimulus = PID_pred(1000:10:end-1000);
d.response = mean2(PID(:,1000:10:end-1000));
p = getModelParameters('hill');
p.A = 1;
p.k = 1;
p.n = 2;
p = FitModel2Data(@hill,d,p);
p_hill = FitModel2Data(@hill,d,p);








