% Convert QP to standard form in MINLPLib, and ensure full rank condition
% In selected MINLPLib instances, all satisfies LB=0
% If you want to handle the case LB!=0, write some additional code for it

function [H, U, D, p, A, b] = standard_form(H, p, A, b, Aeq, beq, LB, UB)

    % Converting to standard form
    m = size(A, 1);
    m_eq = size(Aeq, 1);
    n = length(LB);
    A = [Aeq, zeros(m_eq, m); A, eye(m)];
    b = [beq; b];
    num_ub_slack = 0;
    for i = 1 : n
        if UB(i) < 1e12
            A = [A, zeros(size(A, 1), 1); zeros(1, size(A, 2)), 1];
            A(end, i) = 1;
            b = [b; UB(i)];
            num_ub_slack = num_ub_slack + 1;
        end
    end
    num_slack = m + num_ub_slack;
    H = -0.5 * H;
    p = -0.5 * p;
    H = [H, zeros(n, num_slack); zeros(num_slack, n), zeros(num_slack)];
    p = [p; zeros(num_slack, 1)];
    H = 0.5 * (H + H');
    [U, D] = eig(H);
    
    % Full rank condition
    
    if rank(A) < size(A, 1)
        if check_feasibility(A, b) == 0
            error('Infeasible or something went wrong!');
        end
        At = A';
        [~, basic_row] = rref(At);
        A = At(:, basic_row)';
        b = b(basic_row);
    end
    
    for i = 1 : length(b)
        if abs(b(i)) > 1e5
            b(i) = 1;
            A(i, :) = A(i, :) / b(i);
        end
    end
end

