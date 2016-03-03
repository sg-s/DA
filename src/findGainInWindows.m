% findGainInWindows.m
% 
% created by Srinivas Gorur-Shandilya at 11:56 , 23 February 2016. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.


function  [gain,gain_err,plot_data] = findGainInWindows(ons,offs,pred,resp)

warning off
gain = NaN*ons;
gain_err = NaN*ons;
plot_data(1).x = [];
plot_data(1).y = [];

for i = 1:length(ons)
	try
		[ff,gof]=fit(pred(ons(i):offs(i)),resp(ons(i):offs(i)),'poly1');
		gain(i) = ff.p1;
		gain_err(i) = gof.rsquare;
		plot_data(i).x = pred(ons(i):offs(i));
		plot_data(i).y = resp(ons(i):offs(i));;
	catch
	end
end
warning on
