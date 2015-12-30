% This object is used to represent physiological data.  The object uses the
% design_matrix() class to for the data set and has additional fields for the
% sampling frequency.
%
% :Usage:
% ::
%
%     physio_data: data class for creating a physiological data object
%
% :Inputs:
%
%   **dat:**
%        M x N numeric matrix containing Observations and Variables
%
%   **varname:**
%        Cell array containing variable names.  Must
%                             match number of column in data matrix
%
%   **samplefreq:**
%        Sampling frequency
%
% -------------------------------------------------------------------------
% Current Methods for physio_data (inherits from design_matrix class too)
% -------------------------------------------------------------------------
%
% calc_rate                 : Calculate Rate of peaks
% downsample                : downsample dataset
% filter                    : Filter data
% peakdetect                : Find peaks in data
% physio_data               : class constructor
% plot                      : Plot data
% save                      : Save object as .mat file
% smooth                    : Apply moving average to data
%
% :Example:
% ::
%
%    pulse = physio_data(data(:,2),{'pulse'},settings.fs);
%
% Also see PhysioData_Tutorial.m in Examples
%
% ..
%     -------------------------------------------------------------------------
%     Author and copyright information:
%     -------------------------------------------------------------------------
%     Copyright (C) 2014  Luke Chang
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
%    -------------------------------------------------------------------------
%
%    NOTES:
%    add ability to import .acq files
%    add ability to quickly check all peaks?
% ..

classdef physio_data < design_matrix
    properties
        % inherits properties from design_matrix
        samplefreq = [];
    end
    
    methods
        function obj = physio_data(dat, varname, samplefreq, varargin)
            class constructor % Initialize instance of comp_model
            
            if(nargin > 2)
                
                %Add Data matrix
                try
                    if(~ismatrix(dat) || ~isnumeric(dat) || iscell(dat))
                        error('Make sure input data is a matrix')
                    end
                    if(length(varname) ~= size(dat,2) || ~iscell(varname));
                        error('Make sure the number of variable names corresponds to number of data columns.')
                    end
                    obj.dat = dat;
                    obj.varname = varname;
                catch err
                    error('Make sure input variable names are in a cell array with length equal to number of data columns and data is a matrix.')
                end
                
                %Add variable names
                obj.varname = varname;
                
                %Add Sampling Frequency
                obj.samplefreq = samplefreq;
                
            elseif(nargin > 1)
                try
                    if(~ismatrix(dat) || ~isnumeric(dat) || iscell(dat))
                        error('Make sure input data is a matrix')
                    end
                    if(length(varname) ~= size(dat,2) || ~iscell(varname));
                        error('Make sure the number of variable names corresponds to number of data columns.')
                    end
                    obj.dat = dat;
                    obj.varname = varname;
                catch err
                    error('Make sure input variable names are in a cell array with length equal to number of data columns and data is a matrix.')
                end
                obj.varname = varname;
            elseif(nargin > 0)
                if(~ismatrix(dat) || ~isnumeric(dat) || iscell(dat))
                    error('Make sure input data is a matrix')
                end
                obj.dat = dat;
            else % if nothing initialize empty object
                return
            end
        end
        
        function obj = save(obj,varargin)
            
            % obj = save(obj)
            %
            % -------------------------------------------------------------------------
            % This function saves physio_data class as a .mat file to fname,
            % or user specified path. fullfile(fpath,obj.model)
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % fpath                 : Specify file name path otherwise uses obj.fname
            %
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % data.save('~/MATLAB/HR_Data.mat')
            %
            % data.save('~/MATLAB')
            %
            % -------------------------------------------------------------------------
            
            if nargin >= 1 %use supplied file name
                if ischar(varargin{1})
                    save(varargin{1}, 'obj')
                end
                
            elseif ~isempty(obj.fname) %use obj.fname
                save(obj.fname, 'obj')
            else
                save([obj.varname{1} '_Data.mat'], 'obj')
            end
        end
        
        function f1 = plot(obj,varargin)
            
            % obj = plot(obj)
            %
            % -------------------------------------------------------------------------
            % This function plots physio_data.dat.
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % fpath                 : Specify file name path otherwise uses obj.fname
            %
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % plot(data)
            %
            % -------------------------------------------------------------------------
            
            figure;
            f1 = plot(obj.dat, 'LineWidth', 2);
            
        end
        
        function obj = filter(obj,varargin)
            
            % obj = filter(obj, varargin)
            %
            % -------------------------------------------------------------------------
            % This function applies a zero-order butterworth filter to
            % obj.dat.
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % 'lowpass'            : Followed by lowpass filter cutoff
            % 'highpass'           : Followed by highpass filter cutoff
            % 'bandpass'           : Followed by vector of low and high
            %                        cutoff frequences (e.g., [.05, 2]
            % 'order'              : Follwed by filter order (default 1)
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % data = data.filter('lowpass', .05)
            %
            % data= data.filter('bandpass', [.05, 2])
            %
            % -------------------------------------------------------------------------
            
            % Defaults
            doBandPass = 1;
            doLowPass = 0;
            doHighPass = 0;
            cutoff = [0.015, 2];
            order = 1;
            
            % Parse input
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    if strcmpi('lowpass',varargin{varg})
                        cutoff = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                        doBandPass = 0;
                        doLowPass = 1;
                        doHighPass = 0;
                    elseif strcmpi('highpass',varargin{varg})
                        cutoff = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                        doBandPass = 0;
                        doLowPass = 0;
                        doHighPass = 1;
                    elseif strcmpi('bandpass',varargin{varg})
                        cutoff = varargin{varg + 1};
                        if length(cutoff) ~= 2
                            error('Please include a valid vector of low and high cutoff frequencies')
                        end
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('order',varargin{varg})
                        order = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                    end
                end
            end
            
            nyq = obj.samplefreq / 2;
            
            if doHighPass
                [b, a]=butter(order, cutoff / nyq, 'high');
                obj.dat = filtfilt(b,a,obj.dat); % zero order
            elseif doLowPass
                [b, a]=butter(order, cutoff / nyq, 'low');
                obj.dat = filtfilt(b,a,obj.dat); % zero order
            elseif doBandPass
                [b, a]=butter(order, cutoff / nyq);
                obj.dat = filtfilt(b,a,obj.dat); % zero order
            end
        end
        
        function obj = peakdetect(obj,varargin)
            
            % obj = peakdetect(obj, varargin)
            %
            % -------------------------------------------------------------------------
            % This function applies matlab's peak detection algorithm to
            % data set and adds a vector of the identified peaks to the end of obj.dat
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % 'minPeakHeight'      : Followed by minimum peak height for
            %                        peaks (default: 1 STD above mean)
            %
            % 'plot'               : Plots data with peaks highlighted
            %
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % data = peakdetect(data)
            % data = peakdetect(data,'MinPeakHeight', .5) %Specify minimum peak threshold
            % data = data.peakdetect('plot')
            % -------------------------------------------------------------------------
            
            % NOTES:
            % Could add variable arguments in to findpeaks()
            
            % Defaults
            doMinPeakHeight = 0;
            doPlot = 0;
            %             minPeakHeight = mean(obj.dat) + std(obj.dat); % 1 std above mean
            
            % Parse input
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    if strcmpi('minPeakHeight',varargin{varg})
                        minPeakHeight = varargin{varg + 1};
                        doMinPeakHeight = 1;
                        varargin{varg} = []; varargin{varg + 1} = [];
                    elseif strcmpi('plot',varargin{varg})
                        doPlot = 1;
                        varargin{varg} = [];
                    end
                end
            end
            
            % Use Matlab's findpeaks algorithm to find peaks
            if ~doMinPeakHeight
                [pks,locs] = findpeaks(obj.dat); %returns the indices of the local peaks with no threshold.
            else
                [pks,locs] = findpeaks(obj.dat, 'MinPeakHeight',minPeakHeight); %returns the indices of the local peaks exceeding 1 std of mean of filtered input.
            end
            
            % Add Vector of peaks to end of obj.dat
            obj.varname = horzcat(obj.varname, 'Peaks');
            obj.dat(:,end + 1) = zeros(size(obj.dat));
            obj.dat(locs, end) = 1;
            
            if doPlot
                figure;
                sample = 1:length(obj.dat);
                plot(sample,obj.dat(:, 1:end - 1),sample(locs), pks,'rv','MarkerFaceColor','r');
                xlabel('sample'); ylabel('pulse')
                title('Peak Detection')
            end
        end
        
        function obj = smooth(obj,varargin)
            
            % obj = smooth(obj, varargin)
            %
            % -------------------------------------------------------------------------
            % This function uses matlab's smoothing algorithm to
            % smooth obj.dat with a moving average of X span.
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % 'span'               : Followed by number of samples to span (default: 5)
            %
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % data = smooth(data)
            % data = smooth(data,'span', 100) %Specify the span of the moving average window
            % -------------------------------------------------------------------------
            
            % NOTES:
            % Could add variable arguments in to findpeaks()
            
            % Defaults
            doSpan = 0;
            
            % Parse input
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    if strcmpi('span',varargin{varg})
                        span = varargin{varg + 1};
                        doSpan = 1;
                        varargin{varg} = []; varargin{varg + 1} = [];
                    end
                end
            end
            
            % Use Matlab's findpeaks algorithm to find peaks
            if ~doSpan
                obj.dat = smooth(obj.dat); %span 5
            else
                obj.dat = smooth(obj.dat, span); % custom span
            end
        end
        
        function obj = calc_rate(obj,varargin)
            
            % obj = calc_rate(obj, varargin)
            %
            % -------------------------------------------------------------------------
            % This function calculates the rate of the Peaks calculated with a moving average
            % window length of n samples.
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % 'WindowLength'       : Followed by number of samples of window length
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % data = calc_rate('lowpass', .05)
            %
            % -------------------------------------------------------------------------
            
            % Defaults
            doWindowLen = 0;
            
            % Parse input
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    if strcmpi('WindowLength',varargin{varg})
                        windowlen = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                        doWindowLen = 1;
                    end
                end
            end
            
            % Check if Peaks is in object
            if ~any(strcmpi(obj.varname,'Peaks'))
                error('The variable ''Peaks'' is not in object.  Please run obj.peaks() before running this function.')
            end
            
            obj.dat(:,end + 1) = nan(size(obj,1),1); % Add rate to end of obj.dat
            obj.varname = horzcat(obj.varname,'Rate');
            halfwindow = floor(windowlen/2) * obj.samplefreq;
            i = halfwindow;
            while i + windowlen * obj.samplefreq <= size(obj,1)
                obj.dat(i,end) = sum(obj.dat((i - halfwindow + 1):(i + halfwindow),strcmpi(obj.varname,'Peaks')) / windowlen) * 60;
                i = i + 1;
            end
        end
        
        function obj = downsample(obj,varargin)
            
            % obj = downsample(obj, varargin)
            %
            % -------------------------------------------------------------------------
            % This function downsamples the data.  Can downsample by a
            % factor of N (selects data every N samples) with 'factor' flag.  Alternatively,
            % can 'average' across samples within a grid of time.  This is
            % useful for averaging over TRs for a more reliable estimate when integrating
            % with imaging data.
            %
            % -------------------------------------------------------------------------
            % OPTIONAL INPUTS:
            % -------------------------------------------------------------------------
            % 'Factor'           : Followed by factor to downsample the data.
            % 'Average'          : Followed by window size to average data
            %                      over (e.g., tr for averaging over each tr)
            % -------------------------------------------------------------------------
            % EXAMPLES:
            % -------------------------------------------------------------------------
            % data = data.downsample('Factor', 3) %downsample by a factor of 3
            %
            % data = data.downsample('Average', 1.3) %average over 1.3s TRs
            %
            % -------------------------------------------------------------------------
            
            % Defaults
            doFactor = 0;
            doAverage = 0;
            
            % Parse input
            for varg = 1:length(varargin)
                if ischar(varargin{varg})
                    if strcmpi('Factor',varargin{varg})
                        f = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                        doFactor = 1;
                    elseif strcmpi('Average',varargin{varg})
                        av = varargin{varg + 1};
                        varargin{varg} = []; varargin{varg + 1} = [];
                        doAverage = 1;
                    end
                end
            end
            
            % Downsample by a factor of 'f'
            if doFactor
                obj.dat = downsample(obj.dat, f);
                obj.samplefreq = obj.samplefreq / f;
            end
            
            % Average over window size of 'av' seconds
            if doAverage
                avsamp = av * obj.samplefreq; % average samples
                nTR = floor(size(obj,1)/avsamp); %number of TRs in Dataset
                i = 1;
                start  = 1;
                stop = start + avsamp;
                while i <= nTR
                    ds_dat(i,:) = nanmean(obj.dat(start: stop,:));
                    i = i + 1;
                    start = start + avsamp;
                    stop = stop + avsamp;
                end
                obj.dat = ds_dat; % Update data with averaged data
                obj.samplefreq = 1 / av; % New sampling frequency
            end            
        end
        
    end %methods
end %classdef

