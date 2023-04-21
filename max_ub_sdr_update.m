% Update the upper bound for the maximum value via SDR.

function [ub, y_KKT, I] = max_ub_sdr_update(H, p, A, b, x_lmax, param)

    H_bar = [H p; p' 0];
    A_bar = cell(size(A, 1), 1);
    for i = 1 : size(A, 1)
        A_bar{i} = [zeros(size(A, 2)) 0.5*A(i, :)'; 0.5*A(i, :) 0];
    end
    A_hat = [A zeros(size(A, 1), 1); zeros(1, size(A, 2)) 0];
    B = b * b';
    B_hat = [B zeros(size(B, 1), 1); zeros(1, size(B, 2)) 0];

    [~, t_star] = gurobilp(-ones(size(A, 2), 1), [], [], A, b, zeros(size(A, 2), 1));
    t_star = -t_star;

    if strcmp(param.SOLVER, 'MOSEK')
        [Y, ub] = mosek_sdp(H_bar, A_bar, b, A_hat, B_hat, t_star, param.DNNR_TOL);
    elseif strcmp(param.SOLVER, 'SDPNALP')
        [Y, ub] = sdpnal_plus_sdp(H_bar, A_bar, b, A_hat, B_hat, x_lmax, t_star, param.DNNR_TOL);
    else
        error('Invalid solver name in param.SOLVER');
    end
    
    [V, D] = eig(Y);
    Y_proj = D(end, end) * V(:, end) * V(:, end)';
    y_0 = Y_proj(1:end-1, end);
    y_1 = gurobiqp(eye(size(y_0, 1)), -y_0, A, b);
    [y_KKT, I] = search_local_max_vertex(H, p, A, b, y_1);
   
end