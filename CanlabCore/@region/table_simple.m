function [results_table, results_table_pos, results_table_neg] = table_simple(r)
% Make a simple table of fMRI results with one row per region
%
% [results_table, results_table_pos, results_table_neg] = table_simple(r)
%
% - This is most useful if you've divided an activation map into blobs (regions)
% separated by defined anatomical parcels. 
% - e.g., see image_vector.subdivide_by_atlas
%
% see also region.table, for autolabeling of regions with an atlas and more
% detailed table output.
%
% Tor Wager, April 2021

n_cols = 140;                       % 140 good for HTML reports
sep_str = repmat('_', 1, n_cols);   % see textwrap

dosep = true;
dolegend = true;
table_legend_text = '';

[results_table_pos, results_table_neg] = deal(table());

if isempty(r) || length(r) == 0
    disp('No regions to display');
    return
end

% Separate subregions with positive and negative values if requested
% -------------------------------------------------------------------------
if dosep
    % separate pos and neg
    [poscl, negcl] = posneg_separate(r);
    
    r = [poscl negcl];
    ispos = [true(1, length(poscl)) false(1, length(negcl))]; % logical for splitting combined cl later
    
    clear poscl negcl
    
    fprintf('\n%s\nPositive Effects\n', sep_str)
    
else
    
    fprintf('\n%s\nTable of all regions\n', sep_str)
end


% Volume in mm^3
% count per cubic mm - voxel count * voxel volume
get_vox_volume_in_mm = @(r) size(r.XYZ, 2) .* prod(abs(diag(r.M(1:3, 1:3))));

Volume = zeros(length(r), 1);
for i = 1:length(r)

    Volume(i, 1) = get_vox_volume_in_mm(r(i));
    
end

results_table = table({r.shorttitle}', 'VariableNames', {'Region'});

results_table.Volume = Volume;

results_table.XYZ = round(cat(1, r.mm_center));
    
results_table = [results_table get_signed_max(r, 'Z', 'maxZ')];  % use function because may be empty, handle if so
    
if dosep
    
    results_table_pos = results_table(ispos, :);
    results_table_neg = results_table(~ispos, :);
    
end

% Legend info

voxvol = prod(abs(diag(r(1).M(1:3, 1:3))));

table_legend_text = {'Note: XYZ values are Montreal Neurologic Institute coordinates in mm.'};
table_legend_text{2} = sprintf('Voxel volume is %3.2fmm^3', voxvol);
table_legend_text{3} = 'maxVal is maximum positive or negatively-valued value in the .Z field. Its meaning varies across analyses. For a t-test, it usually reflects the maximum t-value.';

    
%
%     % Sort, if asked for (default = yes)
%     if dosortrows
%     
%         % Replace empty strings so sort will work
%         whempty = cellfun(@isempty, results_table_pos.modal_label_descriptions);
%         results_table_pos.modal_label_descriptions(whempty) = {'X_No_label'};
%         
%         whempty = cellfun(@isempty, results_table_neg.modal_label_descriptions);
%         results_table_neg.modal_label_descriptions(whempty) = {'X_No_label'};
%         
%         results_table_pos = sortrows(results_table_pos, 'modal_label_descriptions');
%         results_table_neg = sortrows(results_table_neg, 'modal_label_descriptions');
% 
%     end
%     
    
% Now split into positive and neg sub-tables and display
if dosep
    % Tables for pos and neg results separately
    
    if any(ispos)
        disp(results_table_pos)
    else
        disp('No regions to display');
    end
    
    fprintf('\nNegative Effects\n')
    if any(~ispos)
        disp(results_table_neg)
    else
        disp('No regions to display');
    end
    
else
    % No separation
    
    if size(results_table, 1) > 0
        disp(results_table)
    else
        disp('No regions to display');
    end
        
end % dosep

if dolegend == false || isempty(table_legend_text)
    return
    
else
    % Print legend text
    
    canlab_print_legend_text(table_legend_text{:});
end

% clean up text for legend
% if length(table_legend_text) > 2
%     table_legend_text(2:3) = [];
% end




end % main function


function val_table = get_signed_max(cl, myfield, tablevarname)
% Returns a table with var "MaxZ" nregions x 1, or empty if cl.Z is empty

%maxZ = @(i) max(double(cl(i).Z));

smax = @(i) signedmax(cl(i).(myfield));

for i = 1:length(cl)
    
    if ~isempty(cl(i).(myfield))
        myZ(i, 1) = smax(i);
    else
        myZ(i, 1) = NaN;
    end
    
    %if isinf(smax(i)), keyboard, end
    
end

% Fix infinite vals - only for .Z . so this is not generalizable beyond
% this function without modifications:
maxZ = norminv(1 - 1E-12);
myZ(myZ > maxZ) = maxZ;

if all(isnan(myZ))
    val_table = [];
else
    val_table = table(myZ,  'VariableNames', {tablevarname});
end

end % function

function val = signedmax(vals)

vals = double(vals);
[maxabs, wh] = max(abs(vals));

val = sign(vals(wh)) .* maxabs;

end

