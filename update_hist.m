% Updating hist info in veri_solve or positive_quadmax

function hist = update_hist(old_hist, time, lb, ub, disp_flag)
    relgap = (ub - lb) / lb;
    hist = [old_hist; time, lb, ub, relgap];
    if disp_flag
        fprintf('%f \t %f \t %f \t %f\n', time, lb, ub, relgap);
    end
end

