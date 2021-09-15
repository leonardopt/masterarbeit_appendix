% Function to print header of the script in the command window 

function layout_print_header(text) %, symbol, rep_symbol, line_break)

%% Set up layout preferences 
% cfg.layout.header.indent      = 20; 
symbol      = '='; 
rep_symbol  = 80;
line_break  = 1;

%% Show header
indent = round((rep_symbol - length(text))/2); 
layout_line(symbol, rep_symbol); layout_line_break(line_break)
layout_indent(indent); fprintf('<strong>%s</strong>\n', text); layout_line_break(line_break);
layout_line(symbol, rep_symbol); layout_line_break(line_break)


end

