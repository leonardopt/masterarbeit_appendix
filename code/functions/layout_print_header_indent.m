% Function to print header of the script in the command window 

function layout_print_header_indent(text, symbol, indent, line_break)

rep_symbol = indent * 2 + length(text);
layout_line(symbol, rep_symbol); layout_line_break(line_break)
layout_indent(indent); fprintf('<strong>%s</strong>\n', text); layout_line_break(line_break);
layout_line(symbol, rep_symbol); layout_line_break(line_break)


end

