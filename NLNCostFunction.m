function [cost] = NLNCostFunction(x,data,CostFunctionHandle)
% unpack data
response = data.response(:);

Rguess = SolveNLNModel(x,data.stimulus);
% calculate the cost
cost = CostFunctionHandle(response(300:end),Rguess(300:end));
