
function [c_sol,x_sol,I,tab_sol]=simplex_two_phases(A,b,c)


m=size(A,1);
n=size(A,2);

%%%%%%%%%%%%%%% normalize data so that the right-hand side vector is nonnegative %%%%%%%%%%%%%%%%%

sgn = sign(b);
pos = find(sgn==0);
sgn(pos) = ones(length(pos), 1);
A=A.*sgn;
b=b.*sgn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('\n\n Phase I.\n\n');

%%%%%%%%%%%%%%%%%%%%%% PHASE I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A2=[A eye(m)];
c2=zeros(n+m,1);
c2(n+1:end)=ones(1,m);
I=n+1:n+m;
n2=n+m;
[tab2,I]=construct_tab(A2,b,c2,n+m,m,I);
max_nb_iters=2000;
debug=0;
[c_sol2,~,I]=simplex_with_Dantzig_pivoting_rule(tab2,I,m,n2,max_nb_iters,debug);
if(max(I)>n)
    if(abs(c_sol2)<1e-10)
        I=I(I<=n);
    else
        error('Problem is infeasible!');
    end
end

[tab_initial,I]=construct_tab(A,b,c,n,m,I);

%%%%%%%%%%%%%%%%%%%%%% END OF PHASE I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%fprintf('\n\n Phase II.\n\n');

%%%%%%%%%%%%%%%%%%%%%% PHASE II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_nb_iters=2000;
debug=0;

[c_sol,x_sol,I,tab_sol,nb_iters]=simplex_with_Dantzig_pivoting_rule(tab_initial,I,m,n,max_nb_iters,debug);

if nb_iters>=max_nb_iters
    %fprintf('Simplex with Dantzig pivoting rule failed, change to Bland pivoting rule');
    [c_sol,x_sol,I,tab_sol,nb_iters]=simplex_with_Bland_pivoting_rule(tab_initial,I,m,n,max_nb_iters*10,debug);
%else
%fprintf('Simplex two Phase succeessfully terminated ');
    if nb_iters>=max_nb_iters*10
        error('Simplex with Bland pivoting rule did not terminate within %d iterations!\n', max_nb_iters*10);
    end
end






end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

