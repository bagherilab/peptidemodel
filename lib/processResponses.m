function s = processResponses(data_file, tab)

% Parse data from excel file.
[num, ~] = xlsread(data_file, tab);

% Field name.
field_name = ['response_' tab];

% Parse numerical data.
s.(field_name).EC50 = num(:,1);
s.(field_name).EC10 = num(:,2);
s.(field_name).UM10 = num(:,3);
s.(field_name).BINARY = num(:,4);

end