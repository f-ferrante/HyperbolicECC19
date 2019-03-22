
%This code solves the optimization problem considered in the paper via a
%bidimensional line search on the scalars mu and alpha. 
%The codes requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
%Any SDP solver can be used.    

clear all;

Lambda =[1 0; 0 2];

H=[0.25 -1;
    0 0.25];

B=eye(2);

N=[1;1];

n=max(size(H));

m=min(size(B));

nd=min(size(N));


mu_min=0.1;
mu_max=5;
alpha_min=0.1;
alpha_max=10;
n_points=40;


mu_v=linspace(mu_min,mu_max,n_points);
alpha_v=linspace(alpha_min,alpha_max,n_points);
gain=[nan,nan];
for i=1:n_points
    
    for j=1:n_points
        
mu=mu_v(i);
alpha=alpha_v(j);

W=diag(sdpvar(n,1));

Gamma=sdpvar(n,n,'symmetric');

Y=sdpvar(m,n,'full');

c=sdpvar(1,1,'full'); 



C1=[-exp(-mu)*W*Lambda, W*H'+Y'*B';
   (W*H'+Y'*B')', -inv(Lambda)*W    
];

C2=[Gamma -N; -N' eye(nd)];

C3=W*(alpha*eye(n)-mu*Lambda)+Gamma;

problem=[C1<=0,C3<=0,C2>=0, W>=1e-9*eye(n),W-c*eye(n)<=0, c>=0];
    
options=sdpsettings('solver','sdpt3','verbose',0);

solution=solvesdp(problem, c,options);

if(solution.problem==0)

W=double(W);

K=double(Y)*inv(W);

c=double(c);
gain(i,j)=max(eig(W))*exp(mu/2);
else
    gain(i,j)=nan;
end

    end
end
%% Seek for the optimal solution
[y,in] = min(gain);

[value,column] = min(y);

[~,row] = min(gain(:,column));

alpha=alpha_v(column);
mu=mu_v(row);

W=diag(sdpvar(n,1));

Gamma=sdpvar(n,n,'symmetric');

Y=sdpvar(m,n,'full');

c=sdpvar(1,1,'full'); 



C1=[-exp(-mu)*W*Lambda, W*H'+Y'*B';
   (W*H'+Y'*B')', -inv(Lambda)*W    
];

C2=[Gamma -N; -N' eye(nd)];

C3=W*(alpha*eye(n)-mu*Lambda)+Gamma;

problem=[C1<=0,C3<=0,C2>=0, W>=1e-9*eye(n),W-c*eye(n)<=0, c>=0];
    
options=sdpsettings('solver','sdpt3','verbose',0);

solution=solvesdp(problem, c,options);

W=double(W);

K=double(Y)*inv(W);

gamma=max(eig(W))*exp(mu/2);

