function [lr_cut, deepened, phi_l_kcut] = lr_cut(D, d, F, f, tuy_cut, konno_cut, bestobjval, phi_0, eta)

    m = size(F, 1);
    l = size(D, 1);
    theta = eta * konno_cut;
    deepened = 0;

    % Constructing matrix for LP relaxation
    % The following code corresponds to the following CVX model:
    % cvx_begin quiet
    %     cvx_solver mosek_2;
    %     variable Lambda_0(l, 1);
    %     variable Lambda(l, m);
    %     variable Lambda_1(l, 1);
    %     variable alph_0;
    %     variable alph(m);
    %     variable b;
    %     minimize(-alph_0 + w' * alph + b);
    %     subject to
    %     Q <= -Lambda_0 * tuy_cut' - tuy_cut * Lambda_0' + Lambda * F + F' * Lambda' + Lambda_1 * theta' + theta * Lambda_1';
    %     2 * d <= 2 * Lambda_0 - 2 * Lambda * w - 2 * Lambda_1 - alph_0 * tuy_cut + F' * alph + b * theta;
    %     Lambda_0 >= 0;
    %     Lambda >= 0;
    %     Lambda_1 >= 0;
    %     alph_0 >= 0;
    %     alph >= 0;
    %     b >= 0;
    % cvx_end
    
    %A = sparse(l * (l + 1), (m + 2) * (l + 1));
    b = zeros(l * (l + 1), 1);
    idx_i = zeros(2 * (m + 2) * l^2 + (m + m + 4) * l, 1);
    idx_j = zeros(2 * (m + 2) * l^2 + (m + m + 4) * l, 1);
    value = zeros(2 * (m + 2) * l^2 + (m + m + 4) * l, 1);
    iter = 1;
    for i = 1 : l
        for j = 1 : l
            b(l * (i-1) + j) = -D(i, j);
            %A(l * (i-1) + j, i) = A(l * (i-1) + j, i) + tuy_cut(j);
            %A(l * (i-1) + j, j) = A(l * (i-1) + j, j) + tuy_cut(i);
            idx_i(iter) = l * (i-1) + j;
            idx_j(iter) = i;
            value(iter) = tuy_cut(j);
            iter = iter + 1;
            idx_i(iter) = l * (i-1) + j;
            idx_j(iter) = j;
            value(iter) = tuy_cut(i);
            iter = iter + 1;
            for k = 1 : m
                %A(l * (i-1) + j, l * k + i) = A(l * (i-1) + j, l * k + i) - F(k, j);
                %A(l * (i-1) + j, l * k + j) = A(l * (i-1) + j, l * k + j) - F(k, i);
                idx_i(iter) = l * (i-1) + j;
                idx_j(iter) = l * k + i;
                value(iter) = - F(k, j);
                iter = iter + 1;
                idx_i(iter) = l * (i-1) + j;
                idx_j(iter) = l * k + j;
                value(iter) = - F(k, i);
                iter = iter + 1;
            end
            %A(l * (i-1) + j, l * (m + 1) + i) = A(l * (i-1) + j, l * (m + 1) + i)-theta(j);
            %A(l * (i-1) + j, l * (m + 1) + j) = A(l * (i-1) + j, l * (m + 1) + j)-theta(i);
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
        %A(l * l + i, i) = -2;
        idx_i(iter) = l * l + i;
        idx_j(iter) = i;
        value(iter) = -2;
        iter = iter + 1;
        for k = 1 : m
            %A(l * l + i, k * l + i) = 2;
            idx_i(iter) = l * l + i;
            idx_j(iter) = k * l + i;
            value(iter) = 2 * f(k);
            iter = iter + 1;
        end
        %A(l * l + i, (m + 1) * l + i) = 2;
        idx_i(iter) = l * l + i;
        idx_j(iter) = (m + 1) * l + i;
        value(iter) = 2;
        iter = iter + 1;
        %A(l * l + i, (m + 2) * l + 1) = tuy_cut(i);
        idx_i(iter) = l * l + i;
        idx_j(iter) = (m + 2) * l + 1;
        value(iter) = tuy_cut(i);
        iter = iter + 1;
        for k = 1 : m
            %A(l * l + i, (m + 2) * l + k + 1) = -F(k, i);
            idx_i(iter) = l * l + i;
            idx_j(iter) = (m + 2) * l + k + 1;
            value(iter) = -F(k, i);
            iter = iter + 1;
        end
        %A(l * l + i, end) = -theta(i);
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