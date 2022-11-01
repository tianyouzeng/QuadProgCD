% Solving the dual of an SDP with the following structure:
% max  <H_bar, Y>;
% s.t. <A_bar(i), Y> = b_i,     i = 1, ..., N;
%      A_hat * Y * A_hat' = B_hat;
%      Y is positively semidefinite;    Y >= 0;    Y(end, end) = 1.
% where A_bar and A_hat has certain relation

function [Y, ub] = mosek_sdp(H_bar, A_bar, b, A_hat, B_hat, t_star, tol)

    k = size(H_bar, 1);
    l = size(A_bar, 1);
    
    % Modeling the SDP
    
    % The scalar part, as in linear optimization examples
    prob.c = [];
    prob.a = sparse([], [], [], length(b) + size(B_hat, 1) + 0.5*k*(k+1) - 1, 0);          % 2 constraints, no scalar variables
    diag_B_hat = diag(B_hat);
    prob.blc = [b; diag_B_hat(1:end-1); zeros(0.5*k*(k+1)-1, 1); 1];                        % Bounds
    prob.buc = [b; diag_B_hat(1:end-1); inf(0.5*k*(k+1)-1, 1); 1];

    % Dimensions of PSD variables
    prob.bardim = k;

    % Coefficients in the objective
    [r1,c1,v1] = find(tril(-H_bar));

    prob.barc.subj = repmat(1,length(v1),1);                % Which PSD variable (j)];           
    prob.barc.subk = r1;                               % Which matrix entry and value ((k,l)->v)
    prob.barc.subl = c1;
    prob.barc.val = v1;

    % Coefficients in the constraints
    prob.bara.subi = [];
    prob.bara.subj = [];
    prob.bara.subk = [];
    prob.bara.subl = [];
    prob.bara.val = [];
    idx = 1;
    for i = 1 : l
        [r1,c1,v1] = find(tril(A_bar{i}));
        prob.bara.subi = [prob.bara.subi; idx * ones(length(v1), 1)];
        prob.bara.subj = [prob.bara.subj; repmat(1,length(v1),1)];
        prob.bara.subk = [prob.bara.subk; r1];
        prob.bara.subl = [prob.bara.subl; c1];
        prob.bara.val = [prob.bara.val; v1];
        idx = idx + 1;
    end
    for i = 1 : l
        [r1,c1,v1] = find(tril(A_hat(i, :)' * A_hat(i, :)));
        prob.bara.subi = [prob.bara.subi; idx * ones(length(v1), 1)];
        prob.bara.subj = [prob.bara.subj; repmat(1,length(v1),1)];
        prob.bara.subk = [prob.bara.subk; r1];
        prob.bara.subl = [prob.bara.subl; c1];
        prob.bara.val = [prob.bara.val; v1];
        idx = idx + 1;
    end
    for i = 1 : k
        for j = 1 : i
            prob.bara.subi = [prob.bara.subi; idx];
            prob.bara.subj = [prob.bara.subj; 1];
            prob.bara.subk = [prob.bara.subk; i];
            prob.bara.subl = [prob.bara.subl; j];
            prob.bara.val = [prob.bara.val; 1];
            idx = idx + 1;
        end
    end

    param.MSK_DPAR_INTPNT_CO_TOL_PFEAS = tol;
    param.MSK_DPAR_INTPNT_CO_TOL_DFEAS = tol;
    param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = tol;
    param.MSK_DPAR_INTPNT_CO_TOL_INFEAS = tol;
    param.MSK_IPAR_AUTO_UPDATE_SOL_INFO = 'MSK_ON';
    %param.MSK_IPAR_INTPNT_SCALING = 'MSK_SCALING_NONE';
    
    [r, res] = mosekopt('minimize echo(0) info', prob, param);

    % Recover the matrix form of primal & dual solution
    Y = zeros(k);
    S = zeros(k);
    lowerTriangleIndices = tril(true(k));
    Y(lowerTriangleIndices) = res.sol.itr.barx; % Fill in the upper piece
    Y = Y + Y.' - diag(diag(Y));
    S(lowerTriangleIndices) = res.sol.itr.bars; % Fill in the upper piece
    S = S + S.' - diag(diag(S));
    gap = abs(res.sol.itr.dobjval - res.sol.itr.pobjval);
    
%     barS = -H_bar;
%     y = res.sol.itr.y;
%     j = 1;
%     for i = 1 : l
%         barS = barS - y(j) * A_bar{i};
%         j = j + 1;
%     end
%     for i = 1 : l
%         barS = barS - y(j) * A_tilde{i};
%         j = j + 1;
%     end
%     R = zeros(k);
%     upperTriangleIndices = triu(true(k));
%     R(upperTriangleIndices) = y(j:end);
%     R = R + R.' - diag(diag(R));
%     barS = barS - R;
    
    % Compute an upper bound from the inexact solution
    % Usually DVIOLCON should be very small
    viol = max([res.info.MSK_DINF_SOL_ITR_DVIOLCON;
        res.info.MSK_DINF_SOL_ITR_PVIOLCON;
        abs(min([min(eig(S)), 0]));
        gap]);
    %viol = 0;% Debug use only, delete this!!!
    ub = (-res.sol.itr.dobjval + (1 + t_star^2) * viol);
end

