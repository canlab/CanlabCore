function canlab_print_legend_text(varargin)
% canlab_print_legend_text(varargin)
% Print a series of strings as a paragraph of text, with delimiters
% before and after the text block. Used in generating HTML reports
% and other functions. Special characters, e.g., line breaks \n, can be
% used.
%
% Example:
% legendtext = sprintf('Wedge plots depict normalized local pattern expression (using the signature weights in the local region), with red indicating positive values and blue negative values. The darker shaded area indicates the standard error of the mean (SEM) across individuals.');
% canlab_print_legend_text('Wedge plot:', legendtext);
%
% Note: If multi-line text is entered and it is split into nonsense letters, try passing 
% in transposed text block

legendtext = varargin; % a series of cells. each cell contains a string

if nargin == 0, return, end  % nothing to display

n_cols = 140;                       % 140 good for HTML reports
sep_str = repmat('_', 1, n_cols);

for i = 1:length(legendtext)
    
    legendtext{i} = sprintf(legendtext{i});  % format line breaks, etc.
    
    legendtext{i} = textwrap(legendtext(i), n_cols);
    legendtext{i} = char(legendtext{i}{:});
    
end

fprintf('\n%s\n', sep_str);

for i = 1:length(legendtext)
    disp(legendtext{i});
end

fprintf('%s\n', sep_str);

end

