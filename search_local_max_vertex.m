% Search for a KKT vertex, starting from a given vertex x_0

function [x_KKT,I] = search_local_max_vertex(H, p, A, b, x_0)
    x = x_0;
    u = 2 * H * x + 2 * p;
    %options = optimoptions('linprog','Display','none', 'Algorithm', 'dual-simplex');
    y = gurobilp(-u, [], [], A, b, zeros(size(x, 1), 1));
    while norm(x - y) > 1e-10 
        x = y;
        u = 2 * H * x + 2 * p;
        y = gurobilp(-u, [], [], A, b, zeros(size(x, 1), 1));
    end
    [~,x_sol,I,~]=simplex_two_phases(A,b,-u);
    % get the position of basic variables for x
    x_KKT = x_sol;
end