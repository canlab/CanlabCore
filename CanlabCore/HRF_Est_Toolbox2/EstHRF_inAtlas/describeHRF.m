function describeHRF(HRF)
    % Generates a description from an HRF Structure.
    % Michael Sun, Ph.D.
    % - Takes HRF Structure object generated from EstimateHRF_inAtlas()
    %
    % *Usage:
    % ::
    %    describeHRF(HRF_structure)

    for c = 1:numel(HRF.params.CondNames)
        for t = 1:numel(HRF.types)
            for r = 1:numel(HRF.region)
            
                disp(['The ', HRF.types{t}, ' fitted ', HRF.region{r}, ' features ', num2str(numel(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phases)), ' phases of activity for the ', HRF.params.CondNames{c}, ' condition']);
                for p = 1:numel(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phases)
                    if HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).peaks>0
                        disp(['Phase ' num2str(p), ', spanning TRs ', num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phases{p}), ' features ', num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).peaks), ' peak(s), (at TR(s) ', num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).time_to_peak), ...
                        ' and features an AUC of ' num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).auc), ' (voxel-normed: ', num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).auc_voxnormed), ').']);
                    end
                    if HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).troughs>0
                     disp(['Phase ' num2str(p), ', spanning TRs ', num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phases{p}), ' features ', num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).troughs), ' trough(s), (at TR(s) ', num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).time_to_peak), ...
                        ' and features an AUC of ' num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).auc), ' (voxel-normed: ', num2str(HRF.fit{t}{r}.(HRF.params.CondNames{c}).phase(p).auc_voxnormed), ').']);
                    end
                end
            
            end
            disp(newline);
        end
    end
end