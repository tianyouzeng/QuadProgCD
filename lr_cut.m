function [lr_cut, deepened, phi_l_kcut] = lr_cut(D, d, F, f, tuy_cut, konno_cut, bestobjval, phi_0, eta)

    m = size(F, 1);
    l = size(D, 1);
    theta = eta * konno_cut;
    deepened = 0;

    % Modeling the LP
    b = zeros(l * (l + 1), 1);
    idx_i = zeros(2 * (m + 2) * l^2 + (m + m + 4) * l, 1);
    idx_j = zeros(2 * (m + 2) * l^2 + (m + m + 4) * l, 1);
    value = zeros(2 * (m + 2) * l^2 + (m + m + 4) * l, 1);
    iter = 1;
    for i = 1 : l
        for j = 1 : l
            b(l * (i-1) + j) = -D(i, j);
            idx_i(iter) = l * (i-1) + j;
            idx_j(iter) = i;
            value(iter) = tuy_cut(j);
            iter = iter + 1;
            idx_i(iter) = l * (i-1) + j;
            idx_j(iter) = j;
            value(iter) = tuy_cut(i);
            iter = iter + 1;
            for k = 1 : m
                idx_i(iter) = l * (i-1) + j;
                idx_j(iter) = l * k + i;
                value(iter) = - F(k, j);
                iter = iter + 1;
                idx_i(iter) = l * (i-1) + j;
                idx_j(iter) = l * k + j;
                value(iter) = - F(k, i);
                iter = iter + 1;
            end
            
            idx_i(iter) = l * (i-1) + j;
            idx_j(iter) = l * (m + 1) + i;
            value(iter) = -theta(j);
            iter = iter + 1;
            idx_i(iter) = l * (i-1) + j;
            idx_j(iter) = l * (m + 1) + j;
            value(iter) = -theta(i);
            iter = iter + 1;
        end
    end
    for i = 1 : l
        b(l * l + i) = - 2 * d(i);
        idx_i(iter) = l * l + i;
        idx_j(iter) = i;
        value(iter) = -2;
        iter = iter + 1;
        for k = 1 : m
            idx_i(iter) = l * l + i;
            idx_j(iter) = k * l + i;
            value(iter) = 2 * f(k);
            iter = iter + 1;
        end
        idx_i(iter) = l * l + i;
        idx_j(iter) = (m + 1) * l + i;
        value(iter) = 2;
        iter = iter + 1;
        idx_i(iter) = l * l + i;
        idx_j(iter) = (m + 2) * l + 1;
        value(iter) = tuy_cut(i);
        iter = iter + 1;
        for k = 1 : m
            idx_i(iter) = l * l + i;
            idx_j(iter) = (m + 2) * l + k + 1;
            value(iter) = -F(k, i);
            iter = iter + 1;
        end
        idx_i(iter) = l * l + i;
        idx_j(iter) = (m + 2) * l + m + 1 + 1;
        value(iter) = -theta(i);
        iter = iter + 1;
    end
    
    A = sparse(idx_i, idx_j, value);

    c = zeros((m + 2) * (l + 1), 1);
    c(end) = 1;
    c((end-m):(end-1)) = f;
    c(end - m - 1) = -1;
    
    % Solving LP relaxation
    [~, fval] = gurobilp(c, A, b, [], [], zeros((m + 2) * (l + 1), 1));

    phi_l_kcut = fval + phi_0;

    if fval <= bestobjval
        lr_cut = theta;
        deepened = 1;
    else
        lr_cut = konno_cut;
        deepened = 0;
    end

end