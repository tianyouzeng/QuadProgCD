% Transform a equality from LP to inequality form

function [D, d, phi_0, F, f, l] = lp_eq2ineq(H, p, A, b, x, I)
    
    n = length(x);

    basic_pos = sort(I);
   
    % positions of nonbasic variables
    non_basic_pos = setdiff(1:size(A, 2), basic_pos);
    non_basic_pos = sort(non_basic_pos);
    
    % computing parameters for coordinate transformation
    % transform the problem into inequality form s.t.
    % vertex x is moved to origin, and
    % the polyhedral is bounded by positive orthant
    B = A(:, basic_pos);
    N = A(:, non_basic_pos);
    F = B \ N;
    f = B \ b;
    H_BB = H(basic_pos, basic_pos);
    H_BN = H(basic_pos, non_basic_pos);
    H_NN = H(non_basic_pos, non_basic_pos);
    p_B = p(basic_pos);
    p_N = p(non_basic_pos);

    D = H_NN + F' * H_BB * F - F' * H_BN - H_BN' * F;
    d = p_N - F' * p_B - F' * H_BB * f + H_BN' * f;
    phi_0 = f' * H_BB * f + 2 * p_B' * f;
    l = size(A, 2) - size(A, 1);
end