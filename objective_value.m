function [v]=objective_value(H,p,x)  

 v=x' * H * x + 2 * p' * x;

end