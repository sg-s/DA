% 
% 
% created by Srinivas Gorur-Shandilya at 10:20 , 09 April 2014. Contact me at http://srinivas.gs/contact/
% 
% This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
% To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
function cost = LCM_Cost_Function(x,d,model)

switch length(fieldnames(d))
	case 0
		error;
	case 1
		error;
	case 2
		error;
	case 3
		K = model(x(1),x(2),x(3));
		Rguess = filter(K,1,d.stimulus);
		Rguess = d.numerator./(1+Rguess);
		cost = Cost2(Rguess(~isnan(Rguess)),d.response(~isnan(Rguess)));
end