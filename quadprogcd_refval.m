% Verifies whether a value c can be attained
% by a convex quadratic function f(x) = x'*H*x + 2*p'*x over the polyhedral set {x : Ax = b, x >= 0}.
% Here, H = U * D * U' for some U, D with positive entries and A, b, p has positive entries.


% Denote by f_* the minimal value of f(x) where x belongs to the
% polyhedral set {x : Ax = b, x >= 0}. Denote by f^* the maximal value of f(x) where x belongs to the
% polyhedral set {x : Ax = b, x >= 0}. 

% Output: status has four possible values: {0,1,2,3}.
% 0 - not sure;
% 1 - attined;
% 2 - not attained;
% 3 - infeasible;

% Output: minval=f_*

% Output: maxval_ub is a value larger than or equal to f^*

% Output: maxval_lb is a value smaller than or equal to f^*

% Requires MOSEK, Gurobi and SDPNAL+ toolbox to run.

function [status, minval, maxval_ub, maxval_lb] = quadprogcd_refval(H, D, U, p, A, b, c, param)
 
    maxtime = param.MAXTIME;
    
    % check if the polyhedral set {x : Ax = b, x >= 0} is empty
    
    if check_feasibility(A, b) == 0
        status = 3;
        minval = Inf;
        maxval_ub = -Inf;
        maxval_lb = -Inf;
        return;
    end

    % compute the minimum value 

    [status, minval] = min_update(H, p, A, b, c);
    if status ~= 0
        maxval_ub = -Inf;
        maxval_lb = -Inf;
        return;
    end
    
    % initializing
    
    lb_global = -Inf;
    ub_global = Inf;
    running_time=0;
    
    tstart = tic;
    if param.PRESOLVE
        ub_global = sdp_lr_presolve(H, p, A, b);
    end

    % performing vertex searching, SDR update and cutting plane methods
    
    while status==0  || running_time>maxtime
        [x_large, lb_1] = search_large_vertex(H, U, p, A, b);  % find a vertex x_large
        lb_global = lb_1;
        if c <= lb_global
            status = 1;
            break;
        end
        
        [x_KKT,I] = search_local_max_vertex(H, p, A, b, x_large);  % find a KKT vertex starting from x_large
        lb_2=objective_value(H,p,x_KKT);
        if(lb_2<lb_1*(1-1e-10))
            error('Something is wrong. The objective value lb_2=%f should be larger than the objective value lb_1=%f. \n\n ', lb_2,lb_1);
        end
        lb_global=max([lb_global; lb_2]);
        
        if c <= lb_global
            status = 1;
            break;
        end
        
        [ub_current, y_KKT, Iy] = max_ub_sdr_update(H, p, A, b, x_KKT, param.SOLVER);
        %if(ub_current>ub_global*(1+1e-10))
        %    error('Something is wrong. The current upper bound ub_current=%f should be smaller than the previous upper bound ub_global=%f.\n\n', ub_current, ub_global);
        %end
        lb_3=objective_value(H,p,y_KKT);
        if lb_3>lb_2
            x_KKT=y_KKT;
            I=Iy;
            lb_global=max([lb_global; lb_3]);
        end
        ub_global = max(lb_global, ub_current);
        

        %%% check if c is attainable %%%
        if c > ub_global
            status = 2;
            break;
        end
        if c <= lb_global
            status = 1;
            break;
        end

        %%%  add konno's cut %%%
        
        % Here, we use Q, d, F, f to store old feasible region
        % We use H, p, A, b to store updated region with Konno's cut
        [Q, d, phi_0, F, f, l] = lp_eq2ineq(H, p, A, b, x_KKT, I);
        HH = H;
        pp = p;
        AA = A;
        bb = b;
        UU = U;
        DD = D;
        [H,U,p,A,b,t_cut,k_cut,status,lb_4] = konno_cut(H,D,U, p, A, b, x_KKT, I, lb_global);
        %[HH,UU,pp,AA,bb,t_cut,k_cut,normineq,status,lb_4] = konno_cut(H,D,U, p, A, b, normineq,x_KKT, I, c);

        %%% check status after adding cut %%%
        if status==1
            lb_global=max(lb_global,lb_4);
        end
        if status==2
            fprintf('\n **After adding Tuy cut, the set becomes infeasible!\n');
              ub_global=c;
        end
        if status==0
            if check_feasibility(A, b) == 0   
                fprintf('\n **After adding Konno cut, the set becomes infeasible!\n');
                ub_global=c;
                status=2;
            end
        end

       %%% add deeper cut %%%

        % We directly pass (H, p, A, b) with Konno's cut to our function
        
        eta = 0.5;
        [l_cut, deepened] = lr_cut(Q, d, F, f, t_cut, k_cut, lb_global - phi_0, phi_0, eta);
        
        if deepened == true
            [H, U, p, A, b] = update_feasible_region(UU, DD, pp, AA, bb, l_cut);
        end
        if deepened == false
            [d_cut, deepened] = dnnr_cut(HH, pp, AA, bb, lb_global, t_cut, k_cut, I, eta);
            if deepened == true
                [H, U, p, A, b] = update_feasible_region(UU, DD, pp, AA, bb, d_cut);
            end
        end
 
        running_time=running_time+toc(tstart);
        
    end
   
    
    maxval_ub = ub_global;
    maxval_lb = lb_global;
end

