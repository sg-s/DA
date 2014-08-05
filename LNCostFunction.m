function cost = LNCostFunction(x,data)
% unpack data
response = data.response(:);

Rguess = SolveLNModel(x,data.stimulus,data.response);

% calculate the cost
cost = Cost2(response(300:end),Rguess(300:end));
