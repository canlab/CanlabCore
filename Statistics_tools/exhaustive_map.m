function [params, maxobj, objmap] = exhaustive_map(p1vals, p2vals, objfun)

    [d1,d2] = meshgrid(p1vals, p2vals);

    objmap = NaN .* zeros(size(d1));

    [m,n] = size(d1);
    for j = 1:m
        for k = 1:n
            objmap(j,k) = objfun(d1(j,k), d2(j,k));
        end
    end
    
    [w1,w2] = find(objmap == max(objmap(:)));
    maxobj = objmap(w1,w2);
    
    for j = 1:size(w1, 1)
        best_d1(j, 1) = d1(w1(j), w2(j));
        best_d2(j, 1) = d2(w1(j), w2(j));
    end
    
    params = [best_d1 best_d2];

    % plot
    
    [data_matrix12, xvals2, yvals2, data_matrix22, bestx, besty, bestz] = surf_plot_tor(objmap, p1vals, p2vals, 'X', 'Y', 'Objective');
end