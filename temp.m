clear;
close all;



g = @circleg;
h_max = 1/4;
T=1;

[p,e,t] = initmesh(g,'hmax',h_max);

tic;
%[xi,u_init, L2E] = RV_solver(p,e,t, T, h_max, @discontinousIC);
[xi,u_init,p,e,tri,t,M] = PDE_Solve_RV(h_max,T, @discontinousIC);
toc;
figure;
pdeplot(p,e,tri,'XYData',u_init, "ZData", u_init);
title("Linear Advection Equation at t = 0 , h_{max} = " + h_max);
xlabel("x");
ylabel("y");

figure;
pdeplot(p,e,tri,'XYData',xi, "ZData", xi);
title("Linear Advection Equation at t = "+ T +" , h_{max} = " + h_max);
xlabel("x");
ylabel("y");

%  convergenceTest(@SUPG_solver, g ,T, [1/4, 1/8, 1/16, 1/32], {@smoothIC, @discontinousIC});
%  convergenceTest(@RV_solver, g ,T, [1/4, 1/8, 1/16, 1/32], {@smoothIC, @discontinousIC});

% plotFigures(@SUPG_solver, g, T, [1/8, 1/16], {@smoothIC, @discontinousIC});
% plotFigures(@RV_solver, g, T, [1/8, 1/16], {@smoothIC, @discontinousIC});

function out = plotFigures(fn, shape, T, h, ICs)
    for i = 1:size(ICs,2)
        for s =1:size(h,2)
            filename = "";
            [p,e,t] = initmesh(shape,'hmax',h(s));
            [xi, u_init, ~] = fn(p,e,t,T,h(s), ICs{i});
            filename = strrep(func2str(fn), "_solver", "") + num2str(h(s)^(-1));
            if func2str(ICs{i}) == "smoothIC"
                filename = filename + "sm";
            elseif func2str(ICs{i}) == "discontinousIC"
                filename = filename + "sh";
            end
            figure;
            pdeplot(p,e,t, 'XYData',u_init, "ZData", u_init);
            title(strrep(func2str(fn),"_","\_") +  " at t = "+ 0 +" , h_{max} = " + h(s) + ", IC = " + func2str(ICs{i}));
            xlabel("x");
            ylabel("y");
            saveas(gcf,"P2/v2/" + filename + "b.png");
            
            figure;
            pdeplot(p,e,t, 'XYData',xi, "ZData", xi);
            title(strrep(func2str(fn),"_","\_") +  " at t = "+ T +" , h_{max} = " + h(s) + ", IC = " + func2str(ICs{i}));
            xlabel("x");
            ylabel("y");
            saveas(gcf,"P2/v2/" + filename + "e.png");
        end
        
    end
end

function y = convergenceTest(fn,shape,T,h, ICs)
    error_vec = zeros(size(ICs,2),size(h,2));
    coeffs = zeros(size(ICs,2), 2);
    figure;
    axes('XScale', 'log', 'YScale', 'log')
    box on
    hold;
    grid;
    loglog(h,h);
    l{1} = "\alpha = 1";
    for i = 1:size(ICs,2)
        for s =1:size(h,2)
            [p,e,t] = initmesh(shape,'hmax',h(s));
            [~, ~, L2E] = fn(p,e,t,T,h(s), ICs{i});
            error_vec(i, s) = L2E;
        end
        coeffs(i,:) = polyfit(log(h), log(error_vec(i,:)),1);
        loglog(h, error_vec(i,:));
        loglog(h, h.^coeffs(i,1).*exp(coeffs(i,2)));
        l{2*i+1} = "Measured for " + func2str(ICs{i});
        l{2*(i+1)} = "\alpha = " + coeffs(i,1) + " for " + func2str(ICs{i});
        
        
    end
    l = {cat(1, l{:})};
    func_name = strrep(func2str(fn), "_", "\_");
    title("Convergence rate of " +func_name);
    xlabel("h_{max} [-]");
    ylabel("L^2-norm error [-]");
    legend(l{1}');
    disp(coeffs);
end

function [xi,u_init, L2E] = SUPG_solver(p,e,t, T, h_max,IC)
    [~,C,M, u_init, beta, ~, area_vec, B, c, bx, by, delta] = assemble(p,e,t, IC);

    delta_C = deltaConvectionAssembler(p,e,t, delta, area_vec, B, c, bx,by);
    delta_Sd = deltaSDAssembler(p,t, delta, area_vec, B, c, bx,by);
    CFL = 0.5;
    k_n = CFL*h_max/max(beta);
    k_n = T/ceil(T/k_n);
    u0 = u_init;
    I = eye(length(p));
    for i = 1:(T/k_n)
        A = (2*(M+delta_C')+k_n*(C+delta_Sd));
        b = (2*(M+delta_C')-k_n*(C+delta_Sd))*u0;
        
        A(e(1,:),:) = I(e(1,:),:);
        b(e(1,:))=0;
        xi = A\b; 
        u0=xi;
        disp(i)
    end
    e = xi - u_init;
    L2E = sqrt(e'*M*e);
end

function [xi, u_init, L2E] = RV_solver(p,e,t, T, h_max, IC)
    [~,C,M, u_init, beta, h, area_vec, B, c, ~, ~, ~] = assemble(p,e,t, IC);
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

        [epsA eps] = epsAassembler(p,e,t, res, xi, area_vec, B, c, beta,h);
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


function [epsA eps_vec] = epsAassembler(p,e,t, res, xi, area_vec, B, c, beta, h)
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
 
        eps_vec(K) = min([Cvel*h(K)*beta_K, Crv*h(K)^2*norm(res(nodes),Inf)/normalization]);
        eps = eps_vec(K);
        
        epsAK = (B(:,K)*B(:,K)' + c(:,K)*c(:,K)')*area_K;
        epsA(nodes,nodes) = epsA(nodes,nodes) +eps*epsAK;
    end
    

end


function delta_C = deltaConvectionAssembler(p,e,t, delta, area_vec, B, c, bx,by)
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

function delta_Sd = deltaSDAssembler(p,t, delta, area_vec, B, c, bx,by)
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

function [A,C,M, u_init, beta, h, area_vec, B, c, bx, by, delta] = assemble(p,e,t,fn)
    N=size(p,2);
    N_t = size(t,2);
    A=sparse(N,N);
    b = zeros(N,1);
    
    C=sparse(N,N);
    M=sparse(N,N);
    u_init = fn(p(1,:),p(2,:))';
    [bx, by] = convectionfield(p(1,:),p(2,:));
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
%         fc = f(xc, yc);
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
%         bK = fc*area_K/3;
        CK = ones(3,1)*(cxmid*B_K+cymid*c_K)'*area_K/3;
        MK = [2 1 1; 1 2 1; 1 1 2]/12*area_K;
       
        A(nodes,nodes) = A(nodes,nodes) +AK;
%         b(nodes)= b(nodes) + bK;
        C(nodes,nodes) = C(nodes,nodes) +CK;
        M(nodes,nodes) = M(nodes,nodes) +MK;
    
        
    
    end
end


function out = smoothIC(x,y)
    r0 = 0.25;
    x0 = 0.3;
    y0 = 0;
    out = 0.5*(1-tanh((((x-x0).^2+(y-y0).^2))./r0^2-1));
end

function out = discontinousIC(x,y)
    r0 = 0.25;
    x0 = 0.3;
    y0 = 0;
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

function [U,U_0,p,e,tri,t,M] = PDE_Solve_RV(hmax,T,u_init)
g = @circleg;
[p,e,tri] = initmesh (g ,'hmax',hmax);
I = eye(length(p));
f_p = f_prim(p(1,:),p(2,:));%convection field
bx = f_p(:,1);
by = f_p(:,2);

%Preallocate matrices
p_a_m = pre_alloc_matrix(p,tri);
M = Mass_assembler_2D(p,tri,p_a_m);
C = ConvectionAssembler2D(p,tri,bx,by,p_a_m);
CFL = 0.5;
k_n = CFL*hmax/norm(f_p,Inf);
[A2,C2,M2, u_init2, beta2, h2, area_vec2, B2, c2, bx2, by2, delta2] = assemble(p,e,tri,@discontinousIC);
U_0 = u_init(p(1,:),p(2,:))';
U_prev = U_0;
assert(isequal(M,M2));
assert(isequal(C,C2));
assert(isequal(U_0, u_init2));
%step once to calculate viscosity
A = (2/k_n)*M+C;
b = ((2/k_n)*M-C)*U_prev;

b2 = ((2/k_n)*M-C)*U_prev;
A2 = (2/k_n)*M+C;

A(e(1,:),:) = I(e(1,:),:);
b(e(1,:)) = 0;

A2(e(1,:),:) = I(e(1,:),:);
b2(e(1,:)) = 0;

U = A\b;
xi = A2\b2;
assert(isequal(U,xi));

t = k_n;
while t < T
    %Want to end at t = T
    %so will at final timestep take a smaller step
    %to reach T exactly
    if t+k_n > T
        k_n = T-t;
    end
    res = M2\(1/k_n*M2*(xi-U_prev) + C*xi);
    
    [epsA eps] = epsAassembler(p,e,tri, res, xi, area_vec2, B2, c2, beta2,h2);
    [Rv, Res] = Rv_Assembler(p,tri,k_n,M,C,U,U_prev,bx,by,p_a_m);
    assert(isequal(res/norm(xi-mean(xi),Inf), Res));
    assert(full(sum(sum(Rv-epsA > eps(1)))));
    A = (2/k_n)*M+C+Rv;
    b = ((2/k_n)*M-C-Rv)*U;
    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:)) = 0;
    U_prev = U;
    U = A\b;
    t = t+k_n;
    disp(t)
end

function f = f_prim(x1,x2)
    f = zeros(length(x1),2);
    for i = 1:length(x1)
        f(i,:) = 2*pi*[-x2(i) x1(i)];
    end
end

function [area,b_Hg,c] = Hat_gradients(x,y)
%HAT_GRADIENTS calculator
area = polyarea(x,y);
b_Hg=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end

function out = pre_alloc_matrix(p,t)
%Returns a matrix consisting of hat gradients
%and area, it is sorted by k
%[area,b,c,h_K]
nt = size(t,2);
out = zeros(nt,8);
for k = 1:nt
    loc2glb = t(1:3,k);
    x = p(1,loc2glb);
    y = p(2,loc2glb);
    [area,b_pam,c] = Hat_gradients(x,y);
    
    %calculated h_K with smallest edge
    edge1 = norm([x(1),y(1)]-[x(2),y(2)]);
    edge2 = norm([x(1),y(1)]-[x(3),y(3)]);
    edge3 = norm([x(2),y(2)]-[x(3),y(3)]);
    h_K = min([edge1,edge2,edge3]);
    
    out(k,:) = [area,b_pam',c',h_K];
end
end

function M = Mass_assembler_2D(p,t,pre_alloc_matrix)
%MASS_ASSEMBLER_2D
%input is point matrix p and connectivity matrix t given by
%initmesh
np = size(p,2);
nt = size(t,2);
M = sparse(np,np);
for k = 1:nt
    loc2glb = t(1:3,k);
    area = pre_alloc_matrix(k,1);
    Mk = [2 1 1;1 2 1;1 1 2]/12*area;
    M(loc2glb,loc2glb) = M(loc2glb,loc2glb)+Mk;
end
end

function C = ConvectionAssembler2D(p,t,bx,by,pre_alloc_matrix)
np=size(p,2);
nt=size(t,2);
C=sparse(np,np);
for k=1:nt
    loc2glb=t(1:3,k);
    area = pre_alloc_matrix(k,1);
    b_C = pre_alloc_matrix(k,2:4)';
    c = pre_alloc_matrix(k,5:7)';
    bxmid=mean(bx(loc2glb));
    bymid=mean(by(loc2glb));
    CK=ones(3,1)*(bxmid*b_C+bymid*c)'*area/3;
    C(loc2glb,loc2glb)=C(loc2glb,loc2glb)+CK;
end
end

function [Rv, Res] = Rv_Assembler(p,t,k_n,M,C,U,U_prev,bx,by,pre_alloc_matrix)
C_vel = 0.25;
C_rv = 1;

%Calculate residual using L2-projection and norm it
Res = M\((1/k_n)*M*(U-U_prev)+C*U);
Res = Res/norm(U-mean(U),Inf);

np = size(p,2);
nt = size(t,2);
Rv = sparse(np,np);
for k = 1:nt
    %Stiffness
    loc2glb = t(1:3,k);
    area = pre_alloc_matrix(k,1);
    b_RV = pre_alloc_matrix(k,2:4)';
    c = pre_alloc_matrix(k,5:7)';
    Ak = (b_RV*b_RV'+c*c')*area;

    %Epsilon
    beta_K = max(sqrt(bx(loc2glb).^2+by(loc2glb).^2));
    h_K = pre_alloc_matrix(k,8);
    eps_K = min([C_vel*h_K*beta_K C_rv*h_K^2*norm(Res(loc2glb),Inf)]);
    
    Rv(loc2glb,loc2glb) = Rv(loc2glb,loc2glb)+eps_K*Ak;
end
end

end

