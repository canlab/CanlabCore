function files = scan_get_files(n, filt, mesg, wd)
% Drop-in replacement for spm_get. 
% Main advantage is that it can use spm_select for SPM5+ to mimic older spm_get behavior.
%
% :Usage:
% ::
%
%     files = scan_get_files(n, filt, mesg, wd)

    if(~exist('wd', 'var') || isempty(wd))
        wd = pwd();
    end

    switch(spm('Ver', [], 1))
        case 'SPM2'
            files = spm_get(n, filt, mesg, wd);
        case {'SPM5', 'SPM8'}
            if(any(n < 0))
                typ = 'dir';
                n = abs(n);
                select_filt = strrep(filt, '*', '.*');
            else
                switch(deblank(filt))
                    case {'*.img' '*img' '*.nii' '*nii' '*.hdr' '*hdr'}
                        typ = 'image';
                        select_filt = '.*';
                    case {'*.mat' '*mat'}
                        typ = 'mat';
                        select_filt = '.*';
                    otherwise
                        typ = 'any';
                        select_filt = strrep(filt, '*', '.*');
                end
            end
            files = spm_select(n, typ, mesg, [], wd, select_filt);
    end
end
