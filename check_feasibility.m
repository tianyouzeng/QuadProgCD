% Feasibility test for the region {Ax = b, x >= 0}

function [is_feasible] = check_feasibility(A, b)
    is_feasible = 1;
    [x, ~, exf] = gurobilp(zeros(size(A, 2), 1), [], [], A, b, zeros(size(A, 2), 1));
    if exf == -2
        is_feasible = 0;
    end
end

