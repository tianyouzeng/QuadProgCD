% Solve the linear relaxation for the SDP stated in mosek_sdp.m, by dropping the semidefinite constraint

function fval = sdp_lr_presolve(H, p, A, b)

    H_bar = [H p; p' 0];
    A_bar = cell(size(A, 1), 1);
    for i = 1 : size(A, 1)
        A_bar{i} = [zeros(size(A, 2)) 0.5*A(i, :)'; 0.5*A(i, :) 0];
    end
    A_hat = [A zeros(size(A, 1), 1); zeros(1, size(A, 2)) 0];
    B = b * b';
    B_hat = [B zeros(size(B, 1), 1); zeros(1, size(B, 2)) 0];
    
    k = size(H_bar, 1);
    l = size(A_bar, 1);
    N = zeros(k, k);
    N(end, end) = 1;
    
    
%     cvx_begin quiet
%         cvx_solver mosek_2;
%         variable Y(k, k) symmetric;
%         maximize trace(H_bar * Y)
%         subject to
%         for i = 1 : l
%             trace(A_bar{i} * Y) == b(i);
%         end
%         %trace(A_h * Y) == b'*b;
%         A_hat * Y * A_hat' == B_hat;
%         %for i = 1 : l
%         %    trace(A_hat(i,:)'*A_hat(i,:)*Y) == B_hat(i, i);
%         %end
%         Y >= 0;
%         trace(N * Y) == 1;
%     cvx_end

    
    f = reshape(H_bar, [], 1);
    Aeq = zeros(2 * l + 1, k*k);
    idx = 1;
    for i = 1 : l
        Aeq(idx, :) = reshape(A_bar{i}, [], 1);
        idx = idx + 1;
    end
    for i = 1 : l
            Aeq(idx, :) = reshape(A_hat(i, :)' * A_hat(i, :), [], 1);
            idx = idx + 1;
    end
    Aeq(end, end) = 1;
    beq = [b; diag(B_hat(1:l, 1:l)); 1];
    [~, fval] = gurobilp(-f, [], [], Aeq, beq, zeros(k*k, 1));
    fval = -fval;
end

