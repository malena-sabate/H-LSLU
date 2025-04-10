function [x_approx,output] = hybrid_LSLU(A,b,x0, x_true, max_iter,RegParam)
% hybrid LSLU codes
% Function to compute Hessenberg form with pivoting for a non-square matrix A
% and use that to solve linear discrete ill posed problems with added
% Tikhonov regualrization on the projected problem.
    
% Initialize values
gstop_vector = zeros(max_iter,1);
d = A'*b;                               % Atransp_times_vec(A, b);
n = length(d);
m = length(b);

% Allocate space
w = zeros(max_iter+1);
h = zeros(max_iter+1);
Enrm =  zeros(max_iter,1);
Rnrm =  zeros(max_iter,1);
RegP =  zeros(max_iter,1);
xk = zeros(n,max_iter);

% Initialize LSLU with Pivoting
r = b-A*x0;
p = (1:m)';     % indeces for solution space
p2 = (1:n)';    % indeces for image space

[~,i0] = max(abs(r));
beta = zeros(max(size(A)),1);
beta(1) = r(i0);
Dk(:,1) = r/beta(1);
% Swap p(1) and p(i0)
p([1 i0]) = p([i0 1]);

for k = 1:max_iter
    q = A'*Dk(:,k);
    for j = 1:k-1
        w(j,k) = q(p2(j));
        q = q - w(j,k) * Lk(:,j);
    end
    if  k < n && sum(q)~=0
        s2 = p2(k:n);
        up2 = q(s2);
        [~,pi02] = max(abs(up2));
        w(k,k) = up2(pi02);
        Lk(:,k) = q/w(k,k);
        %Swaps p2(k) and p2(i02)
        i02 = pi02 + k - 1;
        p2([k i02]) = p2([i02 k]);
    else 
        % Prepare outputs and break
        output.Enrm=Enrm(1:k);
        output.lambda=lambda;
        output.lambda_vector=RegP(1:k);
        output.iteration =k;
        output.Rnrm=Rnrm(1:k);
        output.gstop_vector=gstop_vector;
        output.Lk=Lk;
        output.h=h;
        output.Dk=Dk;
        output.xk=xk;
        break
    end

    u = A * Lk(:,k);
    for j = 1:k
        h(j,k) = u(p(j));
        u = u - h(j,k)* Dk(:,j);
    end
    if  k < m && sum(u)~=0
        s = p(k+1:m);
        up = u(s);
        [~,pi0] = max(abs(up));
        h(k+1,k) = up(pi0);
        Dk(:,k+1) = u/h(k+1,k);

        %Swaps p(k+1) and p(i0)
        i0 = pi0 + k;
        p([k+1 i0]) = p([i0 k+1]);
    else
        % Prepare outputs and break
        output.Enrm=Enrm(1:k);
        output.lambda=lambda;
        output.lambda_vector=RegP(1:k);
        output.iteration =k;
        output.Rnrm=Rnrm(1:k);
        output.gstop_vector=gstop_vector;
        output.Lk=Lk;
        output.h=h;
        output.Dk=Dk;
        output.xk=xk;
        break
    end

    if isscalar(RegParam)
        lambda = RegParam;
        RegP(k) = lambda;
        % Run stopping criterion
        [U,S,~] = svd(h(1:k+1,1:k));
        if size(S,2) > 1
            S = diag(S);
        else
            S = S(1);
        end
        bhat = U' * beta(1:k+1);
        gstop = gcv_stop_lslu(lambda,S,bhat,beta,m,m);
        gstop_vector(k) = gstop;

    elseif strcmp(RegParam,'optimal_tsvd')
        answer = TSVD_Optimal_lslu(h(1:k+1,1:k),beta,x_true, Lk);% h is k+1 x k

    elseif strcmp(RegParam,'optimal_tik')
        opt = @(l) optimal_modified_lslu(l, x_true, x0, h(1:k+1,1:k),beta,Lk,k); % h is k+1 x k
        SS = svds(sparse(h(1:k+1,1:k)),1);% h is k+1 x k
        lambda = fmincon(opt,0,[],[],[],[],0,SS);
        RegP(k) = lambda;
        gstop_vector(k) = 0;

    elseif strcmp(RegParam,'wgcv')
        [U,S,~] = svd(h(1:k+1,1:k)); % h is k+1 x k
        if size(S,2) > 1
            S = diag(S);
        else
            S = S(1);
        end
        bhat = U' * beta(1:k+1);
        mt = size(A,1);
        omega = (k+1)/mt;
        g = @(l)gcv_lslu(l,S,beta,k,bhat,omega); %performs weighted gcv
        lambda = fmincon(g,0,[],[],[],[],0,S(1));
        gstop = gcv_stop_lslu(lambda,S,bhat,beta,mt,size(A,2));
        RegP(k) = lambda;
        gstop_vector(k) = gstop;
    else
        error('Regularization Parameter not valid')
    end

    %Solves the Least Squares Problem
    if strcmp(RegParam,'optimal_tsvd')
        xk(:,k) = x0 + answer;
    else
        y1 = [h(1:k+1,1:k); lambda*eye(k)]\[beta(1:k+1);zeros(k,1)]; % h is k+1 x k
        xk(:,k) = x0 + (Lk(:,1:k) * y1);
    end
    Rnrm(k) = norm(b-A*xk(:,k))/norm(b);
    Enrm(k) = norm(x_true - xk(:,k))/ norm(x_true);
    % compute implicit quasi-residual norm
    % qRnorm(k) = norm(beta(1:k+1) - h*y1)/norm(b);
    % compute augemnted system residual norm
    % ARnrm(k) = norm(b-A*xk(:,k))^2 + lambda^2 * norm(xk(:,k))^2;
end


% Test: check Hessenberg relationship
%   c = A * Lk(:,1:k);
%   d =  Dk * h;
%   e = A'*Dk;
%   f = Lk*w;
% output.c=c;
% output.d=d;
% output.e=e;
% output.f=f;

x_approx = xk(:,k);

% Prepare outputs
output.Enrm=Enrm(1:k);
output.RegP=RegP(1:k);
output.iteration =k;
output.Rnrm=Rnrm(1:k);
output.gstop_vector=gstop_vector;
output.Lk=Lk;
output.h=h;
output.Dk=Dk;
output.xk=xk;
end

%%%% Extra functions needed to set regularization parameters

function [g,S] = gcv_lslu(lambda,S,beta,k,bhat,omega)
%GCV function for Hybrid LSLU - performs weighted gcv or regular gcv
%depending on omega (omega = 1 gives regular gcv; select omega between 0<omega<=1 for weighted)
beta2 = norm(beta)^2;
m = length(bhat);
kbeta = k*beta2;
t1_den = 1./ (S.^2 + lambda^2);
t1 = lambda^2.*t1_den;
t2 = abs(bhat(1:k).*t1).^2;
t3 = sum(abs(bhat(k+1:m)).^2);
t4_num = ((1-omega)*S.^2)+lambda^2;
t4 = t4_num.*t1_den;
g = (kbeta*(sum(t2)+t3))/(1+sum(t4))^2;
end

function gstop = gcv_stop_lslu(Regparamk,S,bhat,beta,m,n)
k = length(S);
beta2 = norm(beta)^2;
nbeta = n*beta2;
t2 = sum(abs(bhat(k+1:end)).^2);
Dk = 1./ (S.^2 + Regparamk^2);
t3 = Regparamk^2 .* Dk ;
t1 = abs(bhat(1:k).*t3).^2;
gstop = (nbeta*(sum(t1) + t2)) / ((m-k) + sum(t3))^2; 
end

function opt = optimal_modified_lslu(Regparamk, x_true, x0,h,beta,PLk,k)
[U,S,V] = svd(h, "econ"); %computes the thin SVD of the given matrix
n = size(h,1);
Sk = diag(S);
bhat = U'* beta(1:n,:);
Dk = Sk.^2 + Regparamk^2;
beta_hat = Sk.* bhat(1:k);
zhat  = beta_hat(1:k)./ Dk;
z = V * zhat;
x = x0 + (PLk(:,1:k)* z);
opt = norm(x-x_true);
end

function out = TSVD_Optimal_lslu(h,beta,x_true,PLk)
[U,S,V] = svd(h, "econ"); %svd of projected matrix h
bhat = U'*beta(1:size(U,1),:); % U' * b
SV = diag(S); %pulls the diagonal entries of singular value matrix
k = size(SV,1);

answer = zeros(size(PLk(:,1)));
e_old = norm(x_true);
for h = 1:k
    answer = answer + PLk(:,1:k)*(bhat(h) / SV(h)) * V(:,h);
    e_new = norm(answer - x_true);
    if e_new < e_old
        out = answer;
    end
    e_old = e_new;
end

end