function [Y, ub] = sdpnal_plus_sdp(H_bar, A_bar, b, A_hat, B_hat, x_lmax, t_star, accuracy_flag)
    k = size(H_bar, 1);
    l = size(A_bar, 1);

    
    % norm_H_bar = norm(H_bar);
    % H_bar = H_bar / norm_H_bar;
    
    print_log = evalc('model = ccp_model();');
    Y = var_sdp(k, k);
    model.add_variable(Y);
    model.minimize(inprod(-H_bar, Y));
    for i = 1 : l
        model.add_affine_constraint(inprod(A_bar{i}, Y) == b(i));
    end
    for i = 1 : l      
        model.add_affine_constraint(inprod(A_hat(i, :)' * A_hat(i, :), Y) == B_hat(i, i));
    end
    model.add_affine_constraint(Y >= 0);
    model.add_affine_constraint(Y(k, k) == 1);
    
    
    model.setparameter('maxiter', 50000, 'printlevel', 2, 'tol', 1e-6, 'stopoption', 1);
    %model.solve();
    
        
    if isempty(x_lmax) % Do not use initial value...
        X_init = [];
    else
        X_init = cell(1);
        X_init{1, 1} = [x_lmax; 1] * [x_lmax; 1]';
    end
  
    
    % Setting initial starting point and solve SDP
    % Codes copied from SDPNAL+
    para = model.info.prob;
    if accuracy_flag == 0
        para.OPTIONS.maxtime = 20;
    end
    % Reshape L, U
    for i = 1:1:para.block
        if strcmp(para.blk{i,1}, 'symm')
            para.blk{i,1} = 'u';
        end
        if para.L{i} == -inf
            para.L{i} = [];
        elseif ~strcmp(para.blk{i,1}, 's')
            para.L{i} = vec(para.L{i});
        end
        if para.U{i} == inf
            para.U{i} = [];
        elseif ~strcmp(para.blk{i,1}, 's')
            para.U{i} = vec(para.U{i});
        end
    end
    print_log = evalc('[obj,X,s,y,S,V,y_bar,v,info_sol,runhist] = sdpnalplus(para.blk,para.At,para.C,para.b,para.L,para.U,para.Bt,para.l,para.u,para.OPTIONS,X_init);');
    
    model.info.issolved = 1;
    for i = 1:1:model.info.prob.block
        if strcmp(model.info.prob.blk{i,1}, 'u') || strcmp(model.info.prob.blk{i,1},'l')
           X{i} = reshape(X{i}, model.info.prob.blkorg{i,2},...
                model.info.prob.blkorg{i,3});
        elseif strcmp(model.info.prob.blk{i,1}, 'symm')
           X{i} = reshape(X{i}, model.info.prob.blkorg{i,2},...
                model.info.prob.blkorg{i,3});
           X{i} = triu(X{i}) + triu(X{i},1)';
        end
    end
    model.info.opt_solution = X;        
    input_data = model.info.prob;
    solution.primal_optimal = X;
    solution.dual_optimal.y = y; 
    solution.dual_optimal.Z1 = S;
    solution.dual_optimal.Z2 = V; 
    model.info.dual_opt = solution.dual_optimal; % new for dual optimal
    solution.info = info_sol;
    %save(model.info.prob.name, 'input_data', 'solution');
    
    Y = X{1,1};
    
    % calculating upper bound using inaccurate solution
    %ub_trX = 2;
    %if size(b_ineq, 1) ~= 0
    %    ub_trX = ub_trX + norm(b_ineq)^2 + norm(A_ineq, 1) * norm(b_ineq, 1) + norm(A_ineq, 1)^2;
    %end
    A_star_y = zeros(size(H_bar));
    for i = 1 : l
        A_star_y = A_star_y + y(i) * A_bar{i};
    end
    for i = 1 : l
        A_star_y = A_star_y + y(i + l) * A_hat(i, :)' * A_hat(i, :);
    end
    N = zeros(k, k);
    N(end, end) = 1;
    A_star_y = A_star_y + y(end) * N;
    S = -H_bar - A_star_y - V{1,1};
    
    pconviol = zeros(2 * l, 1);
    for i = 1 : l
        pconviol(i) = trace(A_bar{i} * Y) - b(i);
    end
    for i = 1 : l      
        pconviol(i + l) = trace(A_hat(i, :)' * A_hat(i, :) * Y) - B_hat(i, i);
    end
    pviol = max(abs(pconviol));
    dviol = abs(min(eig(S)));
    gap = abs(obj(1) - obj(2));
    
    viol = max([pviol; dviol; gap]);

    %ub = (-obj(2) + abs(obj(1) - obj(2)) * (1 + t_star)) * norm_H_bar;
    tracemax = 1 + t_star^2;
    % ub = - (obj(2) - viol * tracemax) * norm_H_bar;
    ub = - (obj(2) - viol * tracemax);

end