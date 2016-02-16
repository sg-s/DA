% findWhiffs.m
% 
% created by Srinivas Gorur-Shandilya at 1:30 , 15 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.

function [whiff_ons,whiff_offs] = findWhiffs(X)

[~,p]=findpeaks(X,'MinPeakProminence',std(X)/5); 
whiff_ons = NaN*p;
whiff_offs = NaN*p;
for i = 1:length(p)
	% cut out a snippet around this peak
	try
		temp = X(p(i)-500:p(i));
		this_whiff_on = p(i) - 500 + find(temp<(X(p(i))/exp(1)),1,'last');
		if isempty(this_whiff_on)
			[~,loc]=min(temp);
			this_whiff_on = p(i) - 500 + loc;
		end
		whiff_ons(i) = this_whiff_on;
		temp = X(p(i):p(i)+500);
		this_whiff_off = p(i) + find(temp<(X(p(i))/exp(1)),1,'first');
		if isempty(this_whiff_off)
			[~,loc]=min(temp);
			this_whiff_off = p(i) - 500 + loc;
		end
		whiff_offs(i) = this_whiff_off;
	catch
	end
end
rm_this = isnan(whiff_offs) | isnan(whiff_ons);
whiff_offs(rm_this) = []; whiff_ons(rm_this)=[];
