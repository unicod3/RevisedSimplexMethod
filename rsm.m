%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function revised solves an LPP using revised simplex method. It uses
% big M method to solve an LPP when thre are >= or = constraints present.
%
% Input:-   1)c : The cost vector or the (row) vector containing co-eficients 
%               of deceision variables in the objective function. It is a
%               required parameter.
%           2)b : The (row) vector containing right hand side constant of
%               the constraints. It is a required parameter.
%           3)a : The coefficient matrix of the left hand side of the
%               constraints. it is a required parameter.
%           4)inq : A (row) vector indicating the type of constraints as 1
%               for >=, 0 for = and -1 for <= constraints. If inq is not
%               supplied then it is by default taken that all constraints
%               are of <= type. It is an optional parameter.
%           5)minimize : This parameter indicates whether the objective
%               function is to be minimized. minimized = 1 indicates
%               a mimization problem and minimization = 0 stands for a
%               maximization problem. By default it is taken as 0. It is an
%               optional parameter.
%
% Example : min     z=-3x1+x2+x3
%           s.t.    x1-2x2+x3<= 11
%               	-4x1+x2+2x3>=3
%                       -2x1+x3=1
%                   x1,x2,x3>=0.
% Solution : c=[-3 1 1];b=[11 3 1];a=[1 -2 1;-4 1 2;-2 0 1];inq=[-1 1 0].
% After supplying these inputs call revised(c,b,a,inq,1).
% 
% This function is able to detect almost all types of
% properties/characteristics present in an LPP such as unbounded solution,
% alternate optima, degenaracy/cycling and infeasibilty. It only fails to
% work when there are redundant constraints present in the problem. However,
% it is rare and can be easily avoided by the user by just
% checking/ensuring that rank(a) should not be less than the number of 
% constraints. As finding rank of big matrices has high complexity, this 
% check has not been given here and it is expected that user would take 
% care of such cases. In such cases usually it is easily seen that some 
% costarints are linearly dependent and hence can be eliminated. Rest of
% the cases show good results. For theory of Revised Simplex method and LPP
% one may see "Numerical Optimization with Applications, Chandra S., Jayadeva,
% Mehra A., Alpha Science Internatinal Ltd, 2009."
%
% This code has been written by Bapi Chatterjee as course assignment in the 
% course Numerical Optimization at Indian Institute of Technology Delhi, 
% New Delhi, India. The author shall be thankful for suggesting any
% error/improvment at bhaskerchatterjee@gmail.com.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function rsm(c,b,a,inq,minimize)
if nargin<3||nargin>5
    fprintf('\nError:Number of input arguments are inappropriate!\n');
else
    n=length(c);m=length(b);j=max(abs(c));
    if nargin<4
        minimize=0;
        inq=-ones(m,1);
    elseif nargin<5
        minimize=0;
    end
    if ~isequal(size(a),[m,n])||m~=length(inq)
        fprintf('\nError: Dimension mismatch!\n');
    else
        if minimize==1
            c=-c;
        end
        count=n;nbv=1:n;bv=zeros(1,m);av=zeros(1,m);
        for i=1:m
            if b(i)<0
                a(i,:)=-a(i,:);
                b(i)=-b(i);
            end
            if inq(i)<0
                count=count+1;
                c(count)=0;
                a(i,count)=1;
                bv(i)=count;
            elseif inq(i)==0
                count=count+1;
                c(count)=-10*j;
                a(i,count)=1;
                bv(i)=count;
                av(i)=count;
            else
                count=count+1;
                c(count)=0;
                a(i,count)=-1;
                nbv=[nbv count];
                count=count+1;
                c(count)=-10*j;
                a(i,count)=1;
                av(i)=count;
                bv(i)=count;
            end
        end
        A=[-c;a];
        B_inv=eye(m+1,m+1);
        B_inv(1,2:m+1)=c(bv);
        x_b=B_inv*[0; b'];
        fprintf('\n.............The initial tablaue................\n')
        fprintf('\t z');disp(bv);
        fprintf('--------------------------------------------------\n')
        disp([B_inv x_b])
        flag=0;count=0;of_curr=0;
        while(flag~=1)
            [s,t]=min(B_inv(1,:)*A(:,nbv));
            y=B_inv*A(:,nbv(t));count=count+1;
            if(any(y(2:m+1)>0))
                fprintf('\n.............The %dth tablaue................\n',count)
                fprintf('\t z');disp(bv);
                fprintf('--------------------------------------------------\n')
                disp([B_inv x_b y])
                if count>1 && of_curr==x_b(1)
                    flag=1;
                    if minimize==1
                        x_b(1)=-x_b(1);
                    end
                    fprintf('\nThe given problem has degeneracy!\n');
                    fprintf('\nThe current objective function value=%d.\n',x_b(1));
                    fprintf('\nThe current solution is:\n');
                    for i=1:n
                        found=0;
                        for j=1:m
                            if bv(j)==i
                                fprintf('x%u = %d\n',i,x_b(1+j));found=1;
                            end
                        end
                        if found==0
                            fprintf('x%u = %d\n',i,0);
                        end
                    end
                else
                    of_curr=x_b(1);
                    if(s>=0)
                        flag=1;
                        for i=1:length(av)
                            for j=1:m
                                if av(i)==bv(j)
                                    fprintf('\nThe given LPP is infeasible!\n');
                                    return
                                end
                            end
                        end
                        if minimize==1
                            x_b(1)=-x_b(1);
                        end
                        fprintf('\nReqiured optimization has been achieved!\n');
                        fprintf('\nThe optimum objective function value=%d.\n',x_b(1));
                        fprintf('\nThe optimum solution is:\n');
                        for i=1:n
                            found=0;
                            for j=1:m
                                if bv(j)==i
                                    fprintf('x%u = %d\n',i,x_b(1+j));found=1;
                                end
                            end
                            if found==0
                                fprintf('x%u = %d\n',i,0);
                            end
                        end
                        if (s==0 && any(y(2:m+1)>0))
                            fprintf('\nThe given problem has alternate optima!\n');
                        end
                    else
                        u=10*j;
                        for i=2:m+1
                            if y(i)>0
                                if (x_b(i)/y(i))<u
                                    u=(x_b(i)/y(i));
                                    v=i-1;
                                end
                            end
                        end
                        temp=bv(v);bv(v)=nbv(t);
                        nbv(t)=temp;
                        E=eye(m+1,m+1);
                        E(:,1+v)=-y/y(1+v);
                        E(1+v,1+v)=1/y(1+v);
                        B_inv=E*B_inv;
                        x_b=B_inv*[0; b'];
                    end
                end
            else
                fprintf('\nThe given problem has unbounded solution\n')
                return
            end
        end
    end
end