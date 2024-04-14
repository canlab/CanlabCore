function tbl = get_regions_at_crosshairs(obj, varargin)
    coords = spm_orthviews('pos')';
    xyz = obj.volInfo.xyzlist;
    xyzmm = (obj.volInfo.mat*[xyz, ones(size(xyz,1),1)]')';
    dist = sum((coords - xyzmm(:,1:3)).^2,2);
    crosshairs_vx_ind = find(dist == min(dist));

    obj = obj.replace_empty();
    
    if ~isempty(obj.probability_maps)
        p = full(obj.probability_maps(crosshairs_vx_ind,:));
        region_ind = find(p>0);
        p = p(region_ind);
    else
        p = 1;
        region_ind = obj.dat(crosshairs_vx_ind);
    end

    if isempty(region_ind)
        tbl = table();
        return
    end

    [labels,labels_2,labels_3,labels_4,labels_5,label_descriptions] = deal(cell(length(region_ind),1));
    for i = 1:length(region_ind)
        labels(i) = obj.labels(region_ind(i));
        if length(obj.labels_2) == length(obj.labels)
            labels_2(i) = obj.labels_2(region_ind(i));
        end
        if length(obj.labels_3) == length(obj.labels)
            labels_3(i) = obj.labels_3(region_ind(i));
        end
        if length(obj.labels_4) == length(obj.labels)
            labels_4(i) = obj.labels_4(region_ind(i));
        end
        if length(obj.labels_5) == length(obj.labels)
            labels_5(i) = obj.labels_5(region_ind(i));
        end
        if length(obj.label_descriptions) == length(obj.labels)
            label_descriptions(i) = obj.label_descriptions(region_ind(i));
        end
    end

    tbl = table(labels, labels_2, labels_3, labels_4, labels_5, label_descriptions, p(:),...
        'VariableNames',{'labels','labels_2','labels_3','labels_4','labels_5','label_descriptions','prob'});

    [~,I] = sort(tbl.prob,'descend');
    tbl = tbl(I,:);

    if ismember(varargin,'display')
        spm_orthviews('caption',1,evalc('disp([tbl.labels, arrayfun(@num2str, tbl.prob, ''UniformOutput'', false)])'))
    end
end