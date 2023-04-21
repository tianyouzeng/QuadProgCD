% Updating hist info in veri_solve or positive_quadmax

function hist = update_hist(old_hist, time, lb, ub, disp_flag, refval)
    relgap = (ub - lb) / lb;
    hist = [old_hist; time, lb, ub, relgap];
    if disp_flag
        if nargin == 6
            fprintf('%f \t %f \t %f \t %f \t %f\n', time, lb, ub, relgap, refval);
        else
            fprintf('%f \t %f \t %f \t %f\n', time, lb, ub, relgap);
        end
    end
end

