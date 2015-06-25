
IGA_data.time = 1e-3*(1:length(mean2(data(i).PID)));
IGA_data.stimulus = mean2(data(i).PID);
IGA_data.prediction = mean2(data(i).LinearFit);
IGA_data.response = mean2(data(i).fA);

% throw out first 5 seconds
IGA_data.time = IGA_data.time(5e3:end);
IGA_data.stimulus = IGA_data.stimulus(5e3:end);
IGA_data.response = IGA_data.response(5e3:end);
IGA_data.prediction = IGA_data.prediction(5e3:end);

% remove trend in stimulus
temp = fit(IGA_data.time(:),IGA_data.stimulus(:),'poly2');
IGA_data.stimulus = IGA_data.stimulus - temp(IGA_data.time) + mean(IGA_data.stimulus);

% fix the gain to be exactly 1
x = IGA_data.prediction;
y = IGA_data.response;
rm_this = isnan(x) | isnan(y);
x(rm_this) = [];
y(rm_this) = [];
temp = fit(x,y,'poly1');
IGA_data.prediction = IGA_data.prediction*temp.p1;

% add the name
IGA_data.name = data(i).original_name;
