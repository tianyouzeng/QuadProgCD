% Search for a vertex with large objective value using the positivity of the problem structure

function [x_large, fval] = search_large_vertex(H, U, p, A, b)
    fval = -Inf;
    for i = 1 : size(U, 2)
        u = U(:, i);
        x = gurobilp(-u, [], [], A, b, zeros(size(u, 1), 1));
        fx = x' * H * x + 2 * p' * x;
        if fx > fval
            fval = fx;
            x_large = x;
        end
    end
end

