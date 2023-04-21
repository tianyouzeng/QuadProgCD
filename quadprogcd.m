% Do the following optimization problem until convergence
% max f(x) = x'*H*x + 2*p'*x over the polyhedral set {x : Ax = b, x >= 0}.
% Here, H = U*D*U' for some U, D with positive entries and A, b, p has positive entries.


% Denote by f_* the minimal value of f(x) where x belongs to the
% polyhedral set {x : Ax = b, x >= 0}. Denote by f^* the maximal value of f(x) where x belongs to the
% polyhedral set {x : Ax = b, x >= 0}. 

% Output: status has four possible values: {0,1,3}.
% 0 - not converged before max_iter reached;
% 1 - converged;
% 3 - infeasible;

% Output: maxval_ub is a value larger than or equal to f^*

% Output: maxval_lb is a value smaller than or equal to f^*

% Requires MOSEK, Gurobi and SDPNAL+ toolbox to run.

function [status, maxval_ub, maxval_lb, hist] = quadprogcd(H, D, U, p, A, b, param)
 
    % Accuracy
    eps = param.RELGAP;
    maxtime = param.MAXTIME;
    hist = [];

    % check if the polyhedral set {x : Ax = b, x >= 0} is empty
    
    if check_feasibility(A, b) == 0
        status = 3;
        maxval_ub = -Inf;
        maxval_lb = -Inf;
        return;
    end

    % initializing
    
    lb_global = -Inf;
    ub_global = Inf;
    running_time=0;
    status = 0;
    
    tstart = tic;
    if param.PRESOLVE
        ub_global = sdp_lr_presolve(H, p, A, b);
    end


    % performing vertex searching, SDR update and cutting plane methods
    
    while status==0  && running_time < maxtime
        
        x_KKT = [];
        lb_2 = -Inf;
        
        [x_large, lb_1] = search_large_vertex(H, U, p, A, b);
        [x_KKT,I] = search_local_max_vertex(H, p, A, b, x_large);
        lb_2=objective_value(H,p,x_KKT);
        if(lb_2<lb_1 - abs(lb_1) * 1e-10)
            error('Something is wrong. The objective value lb_2=%f should be larger than the objective value lb_1=%f. \n\n ', lb_2,lb_1);
        end
        lb_global=max([lb_global; lb_2]);
        
        if param.PRESOLVE
            running_time = toc(tstart);
            hist = update_hist(hist, running_time, lb_global, ub_global, param.DISPLAY);
        end
        
        [ub_current, y_KKT, Iy] = max_ub_sdr_update(H, p, A, b, x_KKT, param);
        lb_3=objective_value(H,p,y_KKT);
        if lb_3>lb_2
            x_KKT=y_KKT;
            I=Iy;
            lb_global=max([lb_global; lb_3]);
        end
        lb_global = max([lb_global; lb_3]);
        
        ub_global = max(lb_global, min(ub_current, ub_global));
        running_time = toc(tstart);
        hist = update_hist(hist, running_time, lb_global, ub_global, param.DISPLAY);

        relgap = (ub_global - lb_global) / lb_global;
        
        %%% check if convergence %%%
        if ub_current <= lb_global || relgap < eps
            status = 1;
            break;
        end

        %%%  add konno's cut %%%

        [Q, d, phi_0, F, f, l] = lp_eq2ineq(H, p, A, b, x_KKT, I);
        HH = H;
        pp = p;
        AA = A;
        bb = b;
        UU = U;
        DD = D;
        [H,U,p,A,b,t_cut,k_cut,cut_status,lb_4] = konno_cut(H,D,U, p, A, b, x_KKT, I, lb_global);

        %%% check status after adding cut %%%
        if cut_status==1
            lb_global=max(lb_global,lb_4);
        end
        if cut_status==2 || (cut_status==0 && check_feasibility(A, b) == 0)
                ub_global = lb_global;
                running_time = toc(tstart);
                hist = update_hist(hist, running_time, lb_global, ub_global, param.DISPLAY);
                status = 1;
                break;
        end

        %%% add deeper cut %%%
        
        if strcmp(param.CUTTYPE, 'ALL')
            eta = 1 / param.CUTDEPTH;
            [l_cut, deepened] = lr_cut(Q, d, F, f, t_cut, k_cut, lb_global - phi_0, phi_0, eta);

            if deepened == true
               [H, U, p, A, b] = update_feasible_region(UU, DD, pp, AA, bb, l_cut, I);
            end
            if deepened == false
                [d_cut, deepened] = dnnr_cut(HH, pp, AA, bb, lb_global, t_cut, k_cut, I, eta, param);
                if deepened == true
                    [H, U, p, A, b] = update_feasible_region(UU, DD, pp, AA, bb, d_cut, I);
                end
            end

            if check_feasibility(A, b) == 0   
                ub_global=lb_global;
                running_time = toc(tstart);
                hist = update_hist(hist, running_time, lb_global, ub_global, param.DISPLAY);
                status = 1;
            end
        end
 
        running_time = toc(tstart);
        hist = update_hist(hist, running_time, lb_global, ub_global, param.DISPLAY);
        
    end
   
    
    maxval_ub = ub_global;
    maxval_lb = lb_global;
end

