% 
% 
% created by Srinivas Gorur-Shandilya at 1:06 , 04 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


%% 
% First, we show that our instatenous gain analysis returns a positive result for synthetic data generated with the DA model, which truly does a "dynamic" gain control mechanism. 

load('Synthetic_data_DA.mat','od')


figure('outerposition',[0 0 1300 400],'PaperUnits','points','PaperSize',[1300 400]); hold on
ax(1) = subplot(1,4,1); hold on
ax(2) = subplot(1,4,2); hold on
ax(3) = subplot(1,4,3); hold on
ax(4) = subplot(1,4,4); hold on

% show the prediction vs. the data
plot(od,ax(1),'ioCurve.firing_rate','data_bin_type','dots','plot_type','mean','showr2',true);

% show the distribution of inst. gain
plot(od,ax(2),'pdf.inst_gain_firing','nbins',100);


plot(od,ax(3:4),'instGainAnalysis.firing_rate.mu','history_lengths',logspace(-2,1,30)*1e3,'history_length',[10 500 1e4],'data_bin_type','dots');


set(ax(3),'XScale','log','YScale','log','XLim',[1e-1 1.5])

suptitle('synthetic Data: DA Model')

prettyFig('fs=13;','FixLogX=true;','lw=1.3;','plw=1.2;');

if being_published
	snapnow
	delete(gcf)
end


%%
% We see that part of the gain control is static ($rho$ does not go to zero at short history times), but there is a distinctive "bump" suggesting that there is gain control at this timescale. 

%% Linear Model
% Now, we do a negative control. We generate synthetic data using a simple linear model. 
