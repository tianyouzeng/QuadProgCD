for i = 1 : 100
    load('C:\Users\Logicwith\OneDrive\Research Projects\Convex Quadratic Maximization\Numerical\Instances\PCQMAX500\pcqmax500_' + string(i) + '.mat')
    S = H ./ (U * D * U');
    scale = S(1, 1);
    D = D * scale;
    save('C:\Users\Logicwith\OneDrive\Research Projects\Convex Quadratic Maximization\Numerical\Instances\PCQMAX500-Normal\pcqmax500_' + string(i) + '.mat', 'H', 'U', 'D', 'p', 'A', 'b', 'c');
end