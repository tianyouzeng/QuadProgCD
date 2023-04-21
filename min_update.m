% Solve the minimum for our quadratic program and refresh the status of our problem

function [status, min] = min_update(H, p, A, b, c)
    options = optimoptions('quadprog','Display','none');
    [~, min, exf] = quadprog(2 * H, 2 * p, [], [], A, b, zeros(1, size(A, 2)), [], [], options);
    
    if exf == -2
        status = 3;
        return;
    end
    if exf ~= 1
        error("min_update: something went wrong Iexitflag)");
    end
    if c < min
        status = 2;
        return;
    end
    status = 0;
end