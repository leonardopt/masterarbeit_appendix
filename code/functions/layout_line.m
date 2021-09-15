% Function to print a chosen symbol in the command window n times 
% Useful to print lines to make code output clearer  

function layout_line(symbol, repetitions)

line = [repmat(symbol, 1, repetitions)];
fprintf('%s\n', line);


end