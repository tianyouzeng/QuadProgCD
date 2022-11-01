function [d_cut, deepened, phi_d_kcut] = dnnr_cut(H, p, A, b, c, tuy_cut, konno_cut, idx, eta, param)

    deepened = 0;

    % Region between Tuy's cut and new cut theta...    
    basic_pos = sort(idx);
    non_basic_pos = setdiff(1:size(A, 2), basic_pos);
    non_basic_pos = sort(non_basic_pos);
    t_cut = zeros(size(A, 2), 1);
    k_cut = zeros(size(A, 2), 1);
    t_cut(non_basic_pos) = tuy_cut;
    k_cut(non_basic_pos) = konno_cut;
    theta = eta * k_cut;
    %t_cut_norm = norm([t_cut;1]);
    %theta_norm = norm([theta;1]);
    A = [A, zeros(size(A, 1), 2); -t_cut', 1, 0; theta', 0, 1];
    b = [b; -1; 1];
    H = [H, zeros(size(H, 1), 2); zeros(2, size(H, 1)), zeros(2, 2)];
    p = [p; 0; 0];

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
    [~, ub] = mosek_sdp(H_bar, A_bar, b, A_hat, B_hat, t_star, param.CUT_TOL);
    %[~, ub] = sdpnal_plus_sdp(H_bar, A_bar, b, A_hat, B_hat, [], t_star, 0);
    
    phi_d_kcut = ub;

    % Acceptance or rejection
    if ub < c
        d_cut = theta;
        deepened = 1;
    else
        d_cut = konno_cut;
        deepened = 0;
    end
end