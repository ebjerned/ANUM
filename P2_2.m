clear;
close all;



g = @circleg;
h_max = 1/16;
T=1;
[p,e,t] = initmesh(g,'hmax',h_max);
[A,C,M,u_init,beta,h] = assemble(p,e,t);
S = SDAssembler(p,t);
CFL = 0.5;
k_n = CFL*h_max/max(beta);
k_n = T/ceil(T/k_n);

u0 = u_init;
I = eye(length(p));

[fx fy] = convectionfield(p(1,:), p(2,:));

figure;
pdeplot(p,e,t,'XYData',u_init, "ZData", u_init);
title("Linear Advection Equation at t = "+ 0*k_n+" , h_{max} = " + h_max);
xlabel("x");
ylabel("y");

% for i= 1:(T/k_n)
% %% RV    
%     if(i == 1)
%         b0 = (M-k_n/2*C)*u0;
%         A0 = (M+k_n/2*C);
%         A(e(1,:),:) = I(e(1,:),:);
%         b(e(1,:))=0;
%         xi = A0\b0;
%         u0=xi;
%     else
%         Cvel = 0.25;
%         Crv = 1;
%         res = M\(1/k_n*M*(xi-u0) + C*xi);
%         
%         normalization = norm(xi-mean(xi), Inf);
%         
%         [Amod eps] = residual_vis(p,e,t,res, normalization);
%         b = (M-k_n/2*C-k_n/2*Amod)*u0;
%         A = (M+k_n/2*C+k_n/2*Amod);
%         A(e(1,:),:) = I(e(1,:),:);
%         b(e(1,:))=0;
%         xi = A\b; 
%   
%         
%         u0=xi;
%     end 
%     
% end
tic;
[xi, L2E] = RV_solver(g, T, h_max);
toc;
figure;
pdeplot(p,e,t,'XYData',xi, "ZData", xi);
title("Linear Advection Equation at t = "+ T +" , h_{max} = " + h_max);
xlabel("x");
ylabel("y");


% e = xi - u_init;
% L2E = sqrt(e'*M*e)
% 
% err_vecSUPG = [0.2645, 0.2726, 0.2781, 0.2637];
% err_vecSUPGdisc = [0.3090, 0.3567, 0.3723, 0.3621];
% err_vecRV = [0.1024, 0.0376, 0.0103, 0.0033];
% err_vecRVdisc = [0.1666, 0.1619, 0.1354,0.1065];
% 
% figure;
% h_steps = [1/4, 1/8, 1/16, 1/32];
% 
% coeff = polyfit(log(h_steps), log(err_vecSUPG),1);
% coeff2 = polyfit(log(h_steps), log(err_vecSUPGdisc),1);
% axes('XScale', 'log', 'YScale', 'log')
% box on
% hold;
% grid;
% loglog(h_steps, h_steps);
% loglog(h_steps, err_vecSUPG);
% loglog(h_steps, err_vecSUPGdisc);
% loglog(h_steps, h_steps.^coeff(1).*exp(coeff(2)));
% loglog(h_steps, h_steps.^coeff2(1).*exp(coeff2(2)));
% legend(["\alpha = 1","SUPG 1.1 Measured", "SUPG 1.3 Measured", "SUPG 1.1 \alpha = "+coeff(1), "SUPG 1.3 \alpha= " + coeff2(1)]);
% title("Convergence rate of different IC for SUPG");
% xlabel("h_{max} [-]");
% ylabel("L^2-norm error [-]");
% 
% figure;
% h_steps = [1/4, 1/8, 1/16, 1/32];
% 
% coeff = polyfit(log(h_steps), log(err_vecRV),1);
% coeff2 = polyfit(log(h_steps), log(err_vecRVdisc),1);
% axes('XScale', 'log', 'YScale', 'log')
% box on
% hold;
% grid;
% loglog(h_steps, h_steps);
% loglog(h_steps, err_vecRV);
% loglog(h_steps, err_vecRVdisc);
% loglog(h_steps, h_steps.^coeff(1).*exp(coeff(2)));
% loglog(h_steps, h_steps.^coeff2(1).*exp(coeff2(2)));
% legend(["\alpha = 1","RV 1.1 Measured", "RV 1.3 Measured", "RV 1.1 \alpha = "+coeff(1), "RV 1.3 \alpha= " + coeff2(1)]);
% title("Convergence rate of different IC for RV");
% xlabel("h_{max} [-]");
% ylabel("L^2-norm error [-]");


function [xi, L2E] = SUPG_solver(shape, T, h_max)
    [p,e,t] = initmesh(shape,'hmax',h_max);
    [A,C,M, u_init, beta, h, area_vec, B, c, bx, by, delta] = assemble(p,e,t);
%     delta = zeros(size(t,2), 1);
    delta_C = deltaConvectionAssembler(p,e,t);
    delta_Sd = deltaSDAssembler(p,t);
    delta_C2 = deltaConvectionAssembler2(p,e,t, delta, area_vec, B, c, bx,by);
    delta_Sd2 = deltaSDAssembler2(p,t, delta, area_vec, B, c, bx,by);
    CFL = 0.5;
    k_n = CFL*h_max/max(beta);
    k_n = T/ceil(T/k_n);
    u0 = u_init;
    I = eye(length(p));
    for i = 1:(T/k_n)
        A = (2*(M+delta_C2')+k_n*(C+delta_Sd2));
        b = (2*(M+delta_C2')-k_n*(C+delta_Sd2))*u0;
        
        A(e(1,:),:) = I(e(1,:),:);
        b(e(1,:))=0;
        xi = A\b; 
        u0=xi;
        disp(i)
    end
    e = xi - u_init;
    L2E = sqrt(e'*M*e);
end

function [xi, L2E] = RV_solver(shape, T, h_max)
    [p,e,t] = initmesh(shape,'hmax',h_max);
    [A,C,M, u_init, beta, h, area_vec, B, c, bx, by, delta] = assemble(p,e,t);
    CFL = 0.5;
    k_n = CFL*h_max/max(beta);
    k_n = T/ceil(T/k_n);
    u0 = u_init;
    b = (M-k_n/2*C)*u0;
    A = (M+k_n/2*C);
    I = eye(length(p));
    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:))=0;
    xi = A\b;
    u0=xi;
    for i = 2:(T/k_n)

        res = M\(1/k_n*M*(xi-u0) + C*xi);

%         [epsA eps] = epsAassembler(p,e,t, res, xi);
        [epsA eps] = epsAassembler2(p,e,t, res, xi, area_vec, B, c, beta);
        b = (M-k_n/2*C-k_n/2*epsA)*u0;
        A = (M+k_n/2*C+k_n/2*epsA);
        
        A(e(1,:),:) = I(e(1,:),:);
        b(e(1,:))=0;
        xi = A\b; 
        u0=xi;
        disp(i)

    end
    e = xi - u_init;
    L2E = sqrt(e'*M*e);

end

function [epsA eps_vec] = epsAassembler(p,e,t, res, xi)
    N = size(p,2);
    epsA = sparse(N,N);
    eps_vec = zeros(N,1);
    Cvel = 0.25;
    Crv = 1.0;
    normalization = norm(xi-mean(xi), Inf);
    
    for K = 1:size(t,2)
       
        nodes = t(1:3, K);
        x=p(1,nodes);
        y=p(2,nodes);
        [bx,by] = convectionfield(x,y);
        
        xc = mean(x);
        yc = mean(y);
        fc = f(xc, yc);

        
        [area_K, B,c] = HatGradients(x,y);
        beta_K = max(sqrt(bx.^2 + by.^2));
        eps_vec(K) = min([Cvel*area_K.*beta_K, Crv*area_K^2*norm(res(nodes),Inf)/normalization]);
        eps = eps_vec(K);
        
        epsAK = (B*B' + c*c')*area_K;
        epsA(nodes,nodes) = epsA(nodes,nodes) +eps*epsAK;
    end
    

end
function [epsA eps_vec] = epsAassembler2(p,e,t, res, xi, area_vec, B, c, beta)
    N = size(p,2);
    epsA = sparse(N,N);
    eps_vec = zeros(N,1);
    Cvel = 0.25;
    Crv = 1.0;
    normalization = norm(xi-mean(xi), Inf);
    
    for K = 1:size(t,2)
       
        nodes = t(1:3, K);        
        area_K = area_vec(K);
        beta_K = beta(K);
 
        eps_vec(K) = min([Cvel*area_K.*beta_K, Crv*area_K^2*norm(res(nodes),Inf)/normalization]);
        eps = eps_vec(K);
        
        epsAK = (B(:,K)*B(:,K)' + c(:,K)*c(:,K)')*area_vec(K);
        epsA(nodes,nodes) = epsA(nodes,nodes) +eps*epsAK;
    end
    

end

function delta_C = deltaConvectionAssembler(p,e,t)
    N=size(p,2);
    
    delta_C=sparse(N,N);
    for K=1:size(t,2)
        nodes=t(1:3,K);
        x=p(1,nodes);
        y=p(2,nodes);
        
        [bx, by] = convectionfield(x, y);
        
        [area, B,c] = HatGradients(x,y);
        h = min(pdist([x',y']));
        delta = 0.5*h/norm([bx,by],2);
        
        cxmid = mean(bx);
        cymid = mean(by);   
        dCK = delta*ones(3,1)*(cxmid*B+cymid*c)'*area/3;
        
        delta_C(nodes,nodes)=delta_C(nodes,nodes)+dCK;
    end

end

function delta_C = deltaConvectionAssembler2(p,e,t, delta, area_vec, B, c, bx,by)
    N=size(p,2);
    
    delta_C=sparse(N,N);
    for K=1:size(t,2)
        nodes=t(1:3,K);         
        cxmid = mean(bx(nodes));
        cymid = mean(by(nodes));   
        dCK = delta(K)*ones(3,1)*(cxmid*B(:,K)+cymid*c(:,K))'*area_vec(K)/3;
        
        delta_C(nodes,nodes)=delta_C(nodes,nodes)+dCK;
    end

end

function [A,C,M, u_init, beta, h, area_vec, B, c, bx, by, delta] = assemble(p,e,t)
    N=size(p,2);
    N_t = size(t,2);
    A=sparse(N,N);
    b = zeros(N,1);
    
    C=sparse(N,N);
    M=sparse(N,N);
    u_init = f(p(1,:),p(2,:))';
    [bx by] = convectionfield(p(1,:),p(2,:));
    beta = zeros(N_t,1);
    area_vec = zeros(N_t,1);
    delta = zeros(N_t,1);
    h = zeros(N_t,1);
    
    B = zeros(3,N_t);
    c = zeros(3,N_t);
    for K = 1:N_t
        %page 243 in book
        nodes = t(1:3, K);
        x=p(1,nodes);
        y=p(2,nodes);
        
        xc = mean(x);
        yc = mean(y);
        fc = f(xc, yc);
        cxmid = mean(bx(nodes));
        cymid = mean(by(nodes));
        
        beta_K = max(sqrt(bx(nodes).^2 + by(nodes).^2));
        beta(K)= beta_K;

        
        [area_K, B_K,c_K] = HatGradients(x,y);
        h(K) = min(pdist([x',y']));
        area_vec(K) = area_K;
        B(:,K) = B_K;
        c(:,K) = c_K;
        delta(K) = 0.5*h(K)/norm([bx(nodes),by(nodes)],2);
        
        AK = (B_K*B_K' + c_K*c_K')*area_K; 
        bK = fc*area_K/3;
        CK = ones(3,1)*(cxmid*B_K+cymid*c_K)'*area_K/3;
        MK = [2 1 1; 1 2 1; 1 1 2]/12*area_K;
       
        A(nodes,nodes) = A(nodes,nodes) +AK;
        b(nodes)= b(nodes) + bK;
        C(nodes,nodes) = C(nodes,nodes) +CK;
        M(nodes,nodes) = M(nodes,nodes) +MK;
    
        
    
    end
    %beta = max(beta_K);
end

function out = f(x,y)
    r0 = 0.25;
    x0 = 0.3;
    y0 = 0;
%       out = 0.5*(1-tanh((((x-x0).^2+(y-y0).^2))./r0^2-1));
    for i = 1:size(x,2)
        if (((x(i)-x0)^2+(y(i)-y0)^2) <= r0^2)
           out(i) = 1;%ones(size(x,2));
        else
            out(i) = 0;%zeros(size(x,2));
        end
    end
    
    
end

function [bx, by] = convectionfield(x,y)
    bx = 2*pi*(-y);
    by = 2*pi*x;
end

function [area,b,c] = HatGradients(x,y)
    area=polyarea(x,y);
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end

function Sd = SDAssembler(p,t)
    np=size(p,2);
    nt=size(t,2);
    Sd=sparse(np,np);
    for i=1:nt
        nodes=t(1:3,i);
        x=p(1,nodes);
        y=p(2,nodes);
        [bx,by] = convectionfield(x,y);
        [area,b,c]=HatGradients(x,y);
        cxmid = mean(bx);
        cymid = mean(by);
        SdK=(cxmid*b+cymid*c)*(cxmid*b+cymid*c)'*area;
        Sd(nodes,nodes)=Sd(nodes,nodes)+SdK;
    end
end

function delta_Sd = deltaSDAssembler(p,t)
    N=size(p,2);
    N_t=size(t,2);
    delta_Sd=sparse(N,N);
    for K=1:N_t
        nodes=t(1:3,K);
        x=p(1,nodes);
        y=p(2,nodes);
        [bx, by] = convectionfield(x,y);
        cxmid = mean(bx);
        cymid = mean(by);
        [area_K, B,c] = HatGradients(x,y);
        h = min(pdist([x',y']));
        delta = 0.5*h/norm([bx,by],2);
        dSdK=delta*(cxmid*B+cymid*c)*(cxmid*B+cymid*c)'*area_K;
        delta_Sd(nodes,nodes)=delta_Sd(nodes,nodes)+dSdK;
    end
end

function delta_Sd = deltaSDAssembler2(p,t, delta, area_vec, B, c, bx,by)
    N=size(p,2);
    N_t=size(t,2);
    delta_Sd=sparse(N,N);
    for K=1:N_t
        nodes=t(1:3,K);
        
        cxmid = mean(bx(nodes));
        cymid = mean(by(nodes));
        
        dSdK=delta(K)*(cxmid*B(:,K)+cymid*c(:,K))*(cxmid*B(:,K)+cymid*c(:,K))'*area_vec(K);
        delta_Sd(nodes,nodes)=delta_Sd(nodes,nodes)+dSdK;
    end
end


