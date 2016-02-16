% ORNData/fitNL.m
% fits a nonlinearity to the data
% 
% created by Srinivas Gorur-Shandilya at 1:34 , 15 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [o] = fitNL(o)

% respect constrains on where to calculate this 
if isempty(o.use_these_trials)
	utt = 1:o.n_trials;
else
	utt = o.use_these_trials;		
end
if isempty(o.use_this_segment)
	uts = 1:length(o.stimulus);
else
	uts = o.use_this_segment;
end	




% check if data is there
if ~isempty(o.firing_rate) && ~isempty(o.firing_projected) && ~isempty(o.stimulus)
	r = o.firing_rate(uts,:);
	p = o.firing_projected(uts,:);

	for i = 1:width(r)
		textbar(i,width(r))
		if ismember(i,utt)
			temp1 = p(:,i); temp2 = r(:,i);
			rm_this = isnan(temp1) | isnan(temp2);
			temp1(rm_this) = []; temp2(rm_this) = [];
			ft = fittype('hillFit(x,A,k,n,x_offset)');
			ff = fit(temp1,temp2,ft,'StartPoint',[max(temp2) mean(temp1) 1 0],'Upper',[1e3 max(temp1) 10 Inf],'Lower',[min(temp2)/2 0 0 -Inf]);
			o.firing_projected(:,i) = ff(o.firing_projected(:,i));
		end
	end

end
if ~isempty(o.LFP) && ~isempty(o.LFP_projected) && ~isempty(o.stimulus)
	r = o.LFP(uts,utt);
	p = o.LFP_projected(uts,utt);

	for i = 1:width(r)
		textbar(i,width(r))
		if ismember(i,utt)
			temp1 = p(:); temp2 = r(:);
			rm_this = isnan(temp1) | isnan(temp2);
			temp1(rm_this) = []; temp2(rm_this) = [];
			ff = fit(temp1(1:100:end),temp2(1:100:end),'rat24','StartPoint',[-2 2 1 -1 -2 .1 1]);
			o.LFP_projected(:,i) = ff(o.LFP_projected(:,i));
		end
	end
end