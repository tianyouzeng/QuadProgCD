% Generate a Kvonno's cut at a given KKT vertex x
% lb: current best fval found on vertices

function [H, U, p, A, b, tuy_cut, konno_cut, status, lb] = konno_cut(H, D0, U,p, A, b, x, I, bestobjval)
    
    n=length(x);
    cut = zeros(n, 1);

    status=0;
    lb=-Inf;

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
    
    
    [~, fval] = gurobilp(-ones(l, 1), F, f, [], [], zeros(l, 1));
    env_splx = -2 * fval;
    

    % generate Tuy's cut
    theta = zeros(l, 1);
    for i = 1 : l
        r = roots([D(i, i) 2*d(i) phi_0-bestobjval]);
        if ~isempty(r)
            theta(i) = max(r);
        else % D(i, i) = d(i) = 0    =>    r can be unbounded
            theta(i) = env_splx;
        end
    end

    tuy_cut = 1./theta;
    tuy_cut(tuy_cut<=1e-10)=0;
    
    [~, ~, exitflag] = gurobilp(zeros(l, 1), [F; 1./theta'], [f; 1], [], [], zeros(l, 1));
    if exitflag == -2 % Tuy's cut lead to a infeasible set
        konno_cut = tuy_cut;
        status=2;
        return;
    end
    
    tau = zeros(l, 1);
    gval_y = zeros(l, 1);
    for i = 1 : l
        y = gurobilp(- (theta(i) * D(i, :)' + d), [F; 1./theta'], [f; 1], [], [], zeros(l, 1));
        gval_y(i) = y' * D * y + 2 * d' * y + phi_0;
    end
    if max(gval_y) > bestobjval
        status = 1;
        konno_cut = tuy_cut;
        lb=max(gval_y);
        return;
    else
        for i = 1 : l
            u = zeros(size(F, 1) + 2, 1);
            u(end) = 1;
            t = ones(l, 1) ./ theta;
            A_ineq = [-F', t, D(i, :)'; f', -1, d(i)];
            b_ineq = [-d; bestobjval - phi_0];
            [~, fval, exf] = gurobilp(-u, A_ineq, b_ineq, [], [], zeros(size(F, 1) + 2, 1));
%             u = [-d; bestobjval - phi_0];
%             t = ones(l, 1) ./ theta;
%             A_ineq = [F -f; -t' 1];
%             b_ineq = zeros(size(F, 1) + 1, 1);
%             A_eq = [D(i, :) d(i)];
%             b_eq = 1;
%             [~, fval, exf] = gurobilp(u, A_ineq, b_ineq, A_eq, b_eq, zeros(l + 1, 1));
            if exf == 1
                tau(i) = -fval;
                %tau(i) = fval;
            else
                if exf == -3 % Problem unbounded
                    tau(i) = env_splx;
                    %tau(i) = 1.1e10;
                else
                    error('something wrong');
                end
            end
        end
    end
    
    cut(non_basic_pos) = 1 ./ tau;
    cut(cut<=1e-10)=0;

    cn=norm([cut;1]);

    A = [A  zeros(size(A, 1), 1); -cut'/cn  1];
    b = [b; -1./cn];
    p = [p; 0];
    U = [U; zeros(1, size(U, 2))];
    H = U * D0 * U';

    konno_cut = 1./tau;
    konno_cut(konno_cut<=1e-10)=0;
    
end