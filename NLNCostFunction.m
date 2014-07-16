function [cost] = NLNCostFunction(x,data)
% unpack data
response = data.response(:);

Rguess = SolveNLNModel(x,data.stimulus);

% calculate the cost
cost = Cost2(response(300:end),Rguess(300:end));
