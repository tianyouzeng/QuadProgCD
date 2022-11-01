% Call Gurobi to solve the following quadratic program:
% max  x'*H*x + 2*p'*x
% s.t. Ax = b, x >= 0

function [x, fval, exitflag] = gurobiqp(H, p, A, b)
    model.Q = sparse(H);
    model.obj = 2 * p;
    model.A = sparse(A);
    model.rhs = b;
    model.sense = '=';
    model.lb = zeros(size(H, 1), 1);
    model.ub = Inf(size(H, 1), 1);
    
    params.OutputFlag = 0;
    
    result = gurobi(model, params);
%     x = result.x;
%     fval = result.objval;
%     if strcmp(result.status, 'OPTIMAL') %#ok<*BDSCA>
%         exf = 1;
%     end
%     if strcmp(result.status, 'INFEASIBLE') || strcmp(result.status, 'INF_OR_UNBD') || strcmp(result.status, 'UNBOUNDED')
%         exf = -2;
%     end
    
% Resolve model if status is INF_OR_UNBD
    if strcmp(result.status,'INF_OR_UNBD')
        params.DualReductions = 0;
        warning('Infeasible or unbounded, resolve without dual reductions to determine...');
        result = gurobi(model,params);
    end

    % Collect results
    x = [];
    output.message = result.status;
    output.constrviolation = [];

    if isfield(result,'x')
        x = result.x;
        if nargout > 3
            slack = model.A*x-model.rhs;
            violA = slack(1:size(A,1));
            violAeq = norm(slack((size(A,1)+1):end),inf);
            viollb = model.lb(:)-x;
            violub = 0;
            if isfield(model,'ub')
                violub = x-model.ub(:);
            end
            output.constrviolation = max([0; violA; violAeq; viollb; violub]);
        end
    end

    fval = [];

    if isfield(result,'objval')
        fval = result.objval;
    end

    if strcmp(result.status,'OPTIMAL')
        exitflag = 1; % converged to a solution
    elseif strcmp(result.status,'UNBOUNDED')
        exitflag = -3; % problem is unbounded
    elseif strcmp(result.status,'ITERATION_LIMIT')
        exitflag = 0; % maximum number of iterations reached
    else
        exitflag = -2; % no feasible point found
    end
end

