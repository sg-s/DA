function [cost] = NLNCostFunction(x,data)
% figure out if K is specified or not
if isfield(data,'K')
	Rguess = SolveNLNModel2(x,data.stimulus(:),data.response(:),data.K);
else
	Rguess = SolveNLNModel2(x,data.stimulus(:),data.response(:));
end

% calculate the cost
cost = Cost2(data.response(300:end),Rguess(300:end));
