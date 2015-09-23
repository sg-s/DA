% EffectOfThreshold.m
% this is an interactive version of the gain analysis function, uses Manipulate.m
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [] = EffectOfThreshold(data)

% fit a LN model
[K,~,filtertime] = FindBestFilter(data.PID,data.ORN,[],'filter_length=199;');
LinearFit = convolve(data.time,data.PID,K,filtertime) + mean(data.ORN);
filtertime = filtertime*3e-3;


% assemble data
s = 300; % when we start for the gain analysis
z = length(data.ORN) - 33; % where we end 
x.response = data.ORN(s:z);
x.LinearFit = LinearFit(s:z);
x.stimulus = data.PID(s:z);
x.time = data.time(s:z);
x.filter_length = 200; 
x.K = K; 
x.filtertime = filtertime;

% make the figure
figure('outerposition',[0 0 1000 900]); hold on
ph = [];
ph(1) = subplot(2,3,1); hold on
set(ph(1),'XLim',[min(filtertime) max(filtertime)])
title(ph(1),'Filter')

ph(2) = subplot(2,3,2:3); hold on
set(ph(2),'XLim',[mean(data.time)-5 mean(data.time)+5])
ph(3) = subplot(2,2,3); hold on; axis square
ph(4) = subplot(2,2,4); hold on;
set(ph(4),'XScale','log')



% use Manipulate to plot everything
x.history_lengths= (3*floor(1000*logspace(-2,1,30)/3))/1e3;
p.t_h = 0.12;
p.frac=  0.33;
p.x0 = 0;
Manipulate('EffectOfThresholdEngine',p,x,[],[],ph);

