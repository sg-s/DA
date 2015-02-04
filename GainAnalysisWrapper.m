% GainAnalysisWrapper.m
% A wrapper for Gain Analysis
% the purpose of this function is to do the gain analysis
% and make all the plots needed with just one line. no more fucking around.
% 
% created by Srinivas Gorur-Shandilya at 4:24 , 04 February 2015. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function [ph] = GainAnalysisWrapper(response,prediction,stimulus,time,example_history_length)

[~,~,~,ct]=FindCorrelationTime(stimulus);
dt = mean(diff(time));
ct = ct*dt;


clear x
x.response = response; 
x.prediction = prediction;
x.stimulus = stimulus; 
x.time = time;
x.filter_length = 299; % what does this even do??
ph = [];

rm_this = [find(isnan(response)) find(isnan(prediction)) ];
x.response(rm_this) = [];
x.prediction(rm_this) = [];
x.stimulus(rm_this) = [];
x.time(rm_this) = [];



history_lengths = (3*floor(1000*logspace(log10(ct),1,30)/3))/1e3;
if nargin < 5
	example_history_length = history_lengths(10);
end

f2=figure('outerposition',[0 0 1000 500],'PaperUnits','points','PaperSize',[1000 500]); hold on
ph(3) = subplot(1,2,1); hold on 
axis square
ph(4) = subplot(1,2,2); hold on

hash = DataHash(x);
cached_data = cache(hash);
if isempty(cached_data)
	[p_LN,l,h] = GainAnalysis4(x,history_lengths,example_history_length,ph);
	cache(hash,p_LN);
	% also cache the example history length
	s=abs(l-h);
	s(p_LN(1,:)>0.05)=NaN;
	[~,loc]=max(s);
	ehl = history_lengths(loc);
	cache(DataHash(p_LN),ehl);

else
	p_LN = cached_data;
	ehl = cache(DataHash(p_LN));
	GainAnalysis4(x,history_lengths,ehl,ph,p_LN);
end

xlabel(ph(3),'Prediction (Hz)')
ylabel(ph(3),'Data (Hz)')
set(ph(4),'XScale','log')
