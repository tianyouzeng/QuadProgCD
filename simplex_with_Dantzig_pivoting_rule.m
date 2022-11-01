function [c_sol,x_sol,I,tab,nb_iters]=simplex_with_Dantzig_pivoting_rule(tab,I,m,n,max_nb_iters,debug)

validate_input(tab,I,m,n); % check if tab is a simplex tableau associated with a basic feasible solution

keep_running=true;
nb_iters=0;



while keep_running&& nb_iters<max_nb_iters
    nb_iters=nb_iters+1;
    if any(tab(1,1:n)>1e-10) %check if there is positive reduced cost coefficient.
        [~,pivot_column] = max(tab(1,1:n)); %yes, find the largest one (Dantzig's pivoting rule)
        if all(tab(2:end,pivot_column)<=0) %check if the corresponding column is negative
            error('problem unbounded. All entries <= 0 in column %d',pivot_column); % yes, the problem is unbounded
        else %otherwise, find pivoting row
            d=tab(2:end,pivot_column);
            f=tab(2:end,end);
            J=find(d>1e-10); %find the indices corresponding to positive entries of d;
            [~,j]=min(f(J)./d(J));
            pivot_row=J(j)+1;
            if debug==1
                fprintf('pivot row is %d\n',pivot_row);
                 fprintf('pivot column is %d\n',pivot_column);
            end
            
            %Apply Gauss-Jordan Exchange
            
            %normalize
            tab(pivot_row,:) = tab(pivot_row,:)/tab(pivot_row,pivot_column);    
            %construct the Q matrix
            d_full=-tab(:,pivot_column);
            d_full(pivot_row)=1;
            Q=eye(m+1,m+1);
            Q(:,pivot_row)=d_full;
            %do elementary row operations
            tab=Q*tab;
            %update the basis index
            I(pivot_row-1)=pivot_column;
           
        end
        if debug==1  %print current basic feasible solution
            fprintf('Simplex Iteration %d \n : current basic feasible solution is\n',nb_iters);
            disp(get_current_x(tab,I,n));
            fprintf('Simplex Iteration %d \n : current objective value is %f.\n',nb_iters,tab(1,end));
        end
    else
        keep_running=false;
        %fprintf('Optimal Solution found in %d simplex iterations.\n',nb_iters);
        if debug==1  %print current basic feasible solution
            fprintf('The optimal basic feasible solution is:\n');
            disp(get_current_x(tab,I,n));
        end
            %fprintf('The optimal objective value found by simplex is %f.\n',tab(1,end));
    end
end
x_sol=get_current_x(tab,I,n);
c_sol=tab(1,end);
end


function current_x = get_current_x(tab, I,n)
current_x = zeros(n,1);
current_x(I)=tab(2:end,end);
end

function validate_input(tab,I,m,n)
if ~ismatrix(tab)
    error('tab must be matrix');
end
if size(tab,1)~=m+1
    error('tab must has %d+1 rows',m);
end
if size(tab,2)~=n+1
    error('tab must has %d+1 columns',n);
end
T=tab(2:end,I);
if max(max(abs(T-eye(m))))>1e-10
    error('The submatrix must be an identity matrix');
end
if min(tab(2:end,:))<-1e-10
    error('the last column of tab must be nonnegative');
end
if (max(abs(tab(1,I)))>1e-10)
    error('the elements in first row and column I of tab must be zero');
end
end