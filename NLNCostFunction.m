function [cost] = NLNCostFunction(x,data)
% unpack data
response = data.response(:);

Rguess = SolveNLNModel2(x,data.stimulus,response);

% calculate the cost
cost = Cost2(response(300:end),Rguess(300:end));
