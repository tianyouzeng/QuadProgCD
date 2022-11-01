
function [tab,I]=construct_tab(A,b,c,n,m,I)
if(length(I)==m)
    J=setdiff(1:n,I);
    B=A(:,I);
    N=A(:,J);
    cB=c(I);
    cN=c(J);
    if(det(B)~=0)
        f=B\b;
        if(min(f)>=-1e-10)
            r=cB'*(B\N)-cN';
            t1=zeros(1,n+1);
            t1(end)=cB'*f;
            t1(J)=r;
            tab=[t1;B\A f];
        else
            fprintf('fmin=%f',min(f));
            error('This basic solution is not feasible!');
        end
    else
        fprintf('det(B)=%f',det(B));
        error('B is not invertible!');
    end
else
    if(rank(A)<m)
        error('A must has rank equal to %d!', m);
    else
        J=setdiff(1:n,I);
        j=1;
        while(length(I)<m && j<=length(J))
            if(rank(A(:,[I j]))>rank(A(:,I)))
                I(end+1)=j;
            end
            j=j+1;
        end
        [tab,I]=construct_tab(A,b,c,n,m,I);
    end
end

end