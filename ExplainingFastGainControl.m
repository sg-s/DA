% ExplainingFastGainControl.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

%% Explaining Fast Gain Control
% We have shown that there is a lot of evidence that ORNs modulate gain as a function of recnetly experienced stimuli. We show that the LN model, despite explaining >95% of the data (typically), fails to account for this fast gain control. Can we develop a model to explain this fast gain control?

%%
% In the following sections, we try fitting a Dynamical Adaptation (DA) model to the data. 


% some parameters
font_size = 20;
marker_size = 10;
marker_size2 = 20;

% load data
load('/local-data/DA-paper/data.mat')
td = 7;

%% Fitting a DA Model to 1-octen-3-ol data
% In this section, we fit a DA model to the example 1-octen-3-ol data: 

disp(data(td).original_name)

% fit model to data
if ~exist('p')
	clear d
	d.stimulus = data(td).PID;
	d.response = data(td).ORN;
	p = FitDAModelToData(d,[3.0114e+04 104 0.03 1.4 2 13 2 -157],[9000 50 0 1e-1 2 1e-2 2 -200],[98000 500 1 20 2 40 2 200]);
	clear d
end

DApred=DA_integrate2(data(td).PID,p);

% fit a LN model (re-use old fits)
load('HowGeneralIsGainAdaptation.mat','Filters','HillFit','filtertime')
LinearFit = mean(data(td).ORN)+convolve(data(td).time,data(td).PID,Filters(:,td),filtertime);
LinearFit(LinearFit<0)=0;
LNpred = hill(HillFit(:,td),LinearFit);

% make a plot showing the data, the LN fit, and the DA fit
fh=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
subplot(2,1,1), hold on
plot(data(td).time,data(td).PID,'k');
set(gca,'XLim',[mean(data(td).time)-2 mean(data(td).time)+2])
ylabel('PID Voltage (V)')
hold off


subplot(2,1,2), hold on
plot(data(td).time,data(td).ORN,'k');
plot(data(td).time,LNpred,'g');
plot(data(td).time,DApred,'r')
set(gca,'XLim',[mean(data(td).time)-2 mean(data(td).time)+2])
ylabel('Firing Rate (Hz)')
xlabel('Time (s)')
legend({'Data','LN Fit','DA Fit'})
PrettyFig;
hold off

snapnow;
delete(fh);




%% Version Info
% The file that generated this document is called:
disp(mfilename)

%%
% and its md5 hash is:
Opt.Input = 'file';
disp(DataHash(strcat(mfilename,'.m'),Opt))