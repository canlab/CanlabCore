function [DAT, behdat] = read_excel(excelfilename, has_header_row)
% Reads each sheet of an Excel input file into a data structure DAT
%
% :Usage:
% ::
%
%     DAT = read_excel(excelfilename, has_header_row [1 or 0])
%
%  - reads multiple sheets
%  - adds fields named with variable names for easy access
%  - uses importdata <matlab internal> to do most of the work
%
% ..
%     Author and copyright information:
%     -------------------------------------------------------------------------
%     Copyright (C) 2015 Tor Wager
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ..
%
% :Input:
%
%   **excelfilename:**
%        Char array with full path of input filename
%
% :Outputs:
%
%   **DAT:**
%        Structure with DAT.(sheetname).(varname) data fields
%        DAT.(sheetname).varnames is cell array of var names
%
%   **behdat:**
%        Data structure as read in by importdata
%
% :Examples:
% ::
%
%    excelfilename = fullfile(basedir, 'data', 'SchulzNewcorn_regressor_demograph.xlsx');
%    has_header_row = 1;
%    [DAT, behdat] = read_excel(excelfilename, has_header_row)
%
% :See also:
% importdata, read_database, read_database2, read_edat_output_2008, read_physio_data
%
% ..
%    Programmers' notes:
%    Created by Tor, 7/2815
%    Not tested extensively yet.
% ..

behdat = importdata(excelfilename); % read the original data

mysheets = fieldnames(behdat.data);

for s = 1:length(mysheets)
% 
% % names
% %-------------------
DAT.(mysheets{s}).varnames = behdat.textdata.(mysheets{s})(1, :);

% % All data
% --------------------


for i = 1:length(DAT.(mysheets{s}).varnames)
    
    fname = genvarname(DAT.(mysheets{s}).varnames{i});
    
    % figure out if is text
    mytext = behdat.textdata.(mysheets{s})((has_header_row+1):end, i); % skip first row if headers.
    istext = ~all(cellfun(@isempty, mytext));
    
    if istext
        
    else  % numeric
        
        DAT.(mysheets{s}).(fname) = behdat.data.(mysheets{s})(:, i);
    end
    
end

end

