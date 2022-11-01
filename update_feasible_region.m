% Update the feasible region by the given cutting planes

function [H, U, p, A, b] = update_feasible_region(U_0, D_0, p_0, A_0, b_0, cut, I)

    % Distinguish whether the input cut is the long cut for original
    % problem or the short cut for canonical problem
    
    if length(cut) == size(p_0, 1)
        long_cut = cut;
    else
        basic_pos = sort(I);

        % positions of nonbasic variables
        non_basic_pos = setdiff(1:length(p_0), basic_pos);
        non_basic_pos = sort(non_basic_pos);

        n = size(cut, 2);
        long_cut = zeros(length(p_0), 1);
        long_cut(non_basic_pos) = cut;
    end
    
    A = [A_0  zeros(size(A_0, 1), 1); -long_cut'  1];
    b = [b_0; -1];
    p = [p_0; 0];
    U = [U_0; zeros(1, size(U_0, 2))];
    H = U * D_0 * U';
end

