clear;
close all;
T = 1;
h = [1/8, 1/16];
ICs = {@shockwaveIC};
field = @convectionfieldNonlinear;
% plotFigures(@RV_RK4_solver, Rectg(-2, -2.5, 2, 1.5) , T, h, ICs, field, 1e-2)
fns = {@SUPG_PIC_solver, @SUPG_RK4_solver,@RV_PIC_solver, @RV_RK4_solver, @GFEM_PIC_solver, @GFEM_RK4_solver};

parfor i = 1:6
    plotFigures(fns{i}, Rectg(-2, -2.5, 2, 1.5), T,h, ICs, field,1e-2);
end



function out = plotFigures(fn, shape, T, h, ICs,field,k_n)
    for i = 1:size(ICs,2)
        for s =1:size(h,2)
            filename = "";
            [p,e,t] = initmesh(shape,'hmax',h(s));
            [xi, u_init, ~] = fn(p,e,t,T,h(s), ICs{i},field,k_n);
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
            light;
            colormap white;
            colorbar off;
            view(40,40)
            saveas(gcf,"P3/" + filename + "b.png");
            
            figure;
            pdeplot(p,e,t, 'XYData',xi, "ZData", xi);
            title(strrep(func2str(fn),"_","\_") +  " at t = "+ T +" , h_{max} = " + h(s) + ", IC = " + func2str(ICs{i}));
            xlabel("x");
            ylabel("y");
            light;
            colormap white;
            colorbar off;
            view(40,40)
            saveas(gcf,"P3/" + filename + "e.png");
        end
        
    end
end

function y = convergenceTest(fn,shape,T,h, ICs, field)
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
            [~, ~, L2E] = fn(p,e,t,T,h(s), ICs{i},field, -1);
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

function [xi,u_init, L2E] = SUPG_PIC_solver(p,e,t, T, h_max,IC, field, k_n)
    [~,C,M, u_init, beta, ~, area_vec, B, c, bx, by, delta] = assemble(p,e,t, IC,field);
    if (k_n == -1)
        CFL = 0.5;
        k_n = CFL*h_max/max(beta);
    end
    k_n = T/ceil(T/k_n); 
    u0 = u_init;

    eps = 1e-1;
    
    xi = u0;
    xi_prev = u0;
    picard_next = u0;
    for i = 1:(T/k_n)
        error = 1;
        counter = 0;
        while ((error > eps) && (counter < 20))
            picard_curr = picard_next;
            [C, delta_C, delta_Sd] = updateSUPGmatrices(p,e,t,delta, area_vec, B, c, picard_curr);
            A = (2/k_n)*(M+delta_C')+C+delta_Sd;
            b = ((2/k_n)*(M+delta_C')-C-delta_Sd)*xi_prev;
            picard_next= A\b;
            picard_next(e(1,:)) = u_init(e(1,:));
            
            error = norm(picard_curr-picard_next);
            disp("Picard error: " + num2str(error))
            counter = counter + 1;
        end
        disp(i)
        xi_prev = xi;
        xi = picard_next;
        
    end
    e = xi - u_init;
    L2E = sqrt(e'*M*e);
end

function [xi, u_init, L2E] = RV_PIC_solver(p,e,t, T, h_max, IC,field,k_n)
    [~,C,M, u_init, beta, h, area_vec, B, c, bx, by, ~] = assemble(p,e,t, IC,field);
    if (k_n == -1)
        CFL = 0.5;
        k_n = CFL*h_max/max(beta);
    end
    k_n = T/ceil(T/k_n); 
    u0 = u_init;
    error = 1;
    eps = 1e-1;
    counter = 0;
    while ((error > eps) && (counter < 20))
        picard_old = u0;
        C = convectionAssemblerNonlinear(p,e,t, area_vec,B, c, picard_old);

        b = (2*M-k_n*C)*u_init;
        A = (2*M+k_n*C);


        
        picard_next = A\b;
        picard_next(e(1,:)) = u_init(e(1,:));
        error = norm(picard_next-picard_old);
        picard_old = picard_next;
        u0 = picard_old;
        disp("Picard Iter, err: " + num2str(error));
        counter = counter +1;
    end
    u0 = u_init;
    error = 1;
    picard_old = u_init;
    picard_curr = picard_next;
    for i = 2:(T/k_n)
        disp(" err: " + num2str(error));
        error = 1;
        while (error > eps)
            C = convectionAssemblerNonlinear(p,e,t, area_vec,B, c, picard_curr);
            res = M\(1/k_n*M*(picard_curr-u0) + C*picard_old);

            [epsA, ~] = epsAassembler(p,e,t, res, picard_curr, area_vec, B, c, beta,h);
            b = (2*M - k_n*(C + epsA))*u0;
            A = (2*M +k_n*(C + epsA));


            picard_old = picard_curr;
            picard_next = A\b; 
            picard_next(e(1,:)) = u_init(e(1,:));
            error = norm(picard_next - picard_curr);
            
            picard_curr = picard_next;
            
            disp("Picard Iter, err: " + num2str(error));
        end
        u0 = picard_old;
        disp(i)

    end
    e = picard_curr - u_init;
    L2E = sqrt(e'*M*e);
    xi = picard_curr;

end

function [xi, u_init, L2E] = GFEM_PIC_solver(p,e,t, T, h_max, IC,field,k_n)
    [~,C,M, u_init, beta, ~, area_vec, B, c, ~, ~, ~] = assemble(p,e,t, IC,field);
    if (k_n == -1)
        CFL = 0.5;
        k_n = CFL*h_max/max(beta)*0.1;
    end
    k_n = T/ceil(T/k_n); 
    u0 = u_init;
    xi_curr = u0;
    for i = 1:(T/k_n)
        error = 1;
        eps = 1e-1;
        picard_curr = xi_curr;
        counter = 0;
        while ((error > eps) && (counter < 20))
            C = convectionAssemblerNonlinear(p,e,t, area_vec,B, c, picard_curr);

            b = (2*M-k_n*C)*xi_curr;
            A = (2*M+k_n*C);
            picard_next = A\b;
            picard_next(e(1,:)) = u_init(e(1,:));
            error = norm(picard_next - picard_curr);
            picard_curr = picard_next;
            disp("Picard Iter, err: " + num2str(error));
            counter = counter +1;
        end
        
        xi_curr = picard_curr;
        disp(i)

    end
    xi = xi_curr;
    e = xi - u_init;
    L2E = sqrt(e'*M*e);


end

function [xi,u_init, L2E] = SUPG_RK4_solver(p,e,t, T, h_max,IC, field,k_n)
    [~,C,M, u_init, beta, ~, area_vec, B, c, bx, by, delta] = assemble(p,e,t, IC,field);
    if (k_n == -1)
        CFL = 0.5;
        k_n = CFL*h_max/max(beta)*0.1;
    end
    k_n = T/ceil(T/k_n); 
    
    u0 = u_init;
    for i = 1:(T/k_n)

        k1 = (M+deltaConvectionAssemblerNonlinear(p,e,t, delta, area_vec, B, c ,u0)')\((-convectionAssemblerNonlinear(p,e,t,area_vec, B,c,u0)-deltaSDAssemblerNonlinear(p,t, delta, area_vec, B, c, u0))*u0);
        k2 = (M+deltaConvectionAssemblerNonlinear(p,e,t, delta, area_vec, B, c, u0+k_n/2*k1)')\((-convectionAssemblerNonlinear(p,e,t,area_vec, B,c,u0+k_n/2*k1)-deltaSDAssemblerNonlinear(p,t, delta, area_vec, B, c, u0+k_n/2*k1))*(u0+k_n/2*k1));
        k3 = (M+deltaConvectionAssemblerNonlinear(p,e,t, delta, area_vec, B, c, u0+k_n/2*k2)')\((-convectionAssemblerNonlinear(p,e,t,area_vec, B,c,u0+k_n/2*k2)-deltaSDAssemblerNonlinear(p,t, delta, area_vec, B, c, u0+k_n/2*k2))*(u0+k_n/2*k2));
        k4 = (M+deltaConvectionAssemblerNonlinear(p,e,t, delta, area_vec, B, c, u0+k_n*k3)')\((-convectionAssemblerNonlinear(p,e,t,area_vec, B,c,u0+k_n*k3)-deltaSDAssemblerNonlinear(p,t, delta, area_vec, B, c, u0+k_n*k3))*(u0+k_n*k3));
        
        xi = u0 + k_n/6*(k1+2*k2+2*k3 + k4);
        xi(e(1,:)) = u_init(e(1,:));
        u0 = xi;
        disp(i)
    end
    e = xi - u_init;
    L2E = sqrt(e'*M*e);
end


function [xi, u_init, L2E] = RV_RK4_solver(p,e,t, T, h_max, IC, field,k_n)
    [~,C,M, u_init, beta, h, area_vec, B, c, bx, by, ~] = assemble(p,e,t, IC,field);
    if (k_n == -1)
        CFL = 0.5;
        k_n = CFL*h_max/max(beta)*0.1;
    end
    k_n = T/ceil(T/k_n); 
   
    u0 = u_init;
    
    k1 = M\(-convectionAssemblerNonlinear(p,e,t, area_vec,B, c, u0)*u0);
    k2 = M\(-convectionAssemblerNonlinear(p,e,t, area_vec,B, c, u0+k_n/2*k1)*(u0+k_n/2*k1));
    k3 = M\(-convectionAssemblerNonlinear(p,e,t, area_vec,B, c, u0+k_n/2*k2)*(u0+k_n/2*k2));
    k4 = M\(-convectionAssemblerNonlinear(p,e,t, area_vec,B, c, u0+k_n*k3)*(u0+k_n*k3));
    xi = u0 + k_n/6*(k1+2*k2+2*k3 + k4);
    xi(e(1,:)) = u_init(e(1,:));
  

    for i = 2:(T/k_n)

        res = M\(1/k_n*M*(xi-u0) + C*xi);
        k1 = M\((-convectionAssemblerNonlinear(p,e,t, area_vec, B, c, xi)-epsAassembler(p,e,t,res,xi,area_vec,B,c,beta,h))*xi);
        res = M\(1/k_n*M*((xi+k_n/2*k1)-u0) + C*(xi+k_n/2*k1));
        k2 = M\((-convectionAssemblerNonlinear(p,e,t, area_vec, B, c, xi + k_n/2*k1)-epsAassembler(p,e,t,res, xi+k_n/2*k1, area_vec, B,c, beta,h))*(xi+k_n/2*k1));
        res = M\(1/k_n*M*((xi+k_n/2*k2)-u0) + C*(xi+k_n/2*k2));
        k3 = M\((-convectionAssemblerNonlinear(p,e,t, area_vec, B, c, xi + k_n/2*k2)-epsAassembler(p,e,t,res, xi+k_n/2*k2, area_vec, B,c, beta,h))*(xi+k_n/2*k2));
        res = M\(1/k_n*M*((xi+k_n*k3)-u0) + C*(xi+k_n*k3));
        k4 = M\((-convectionAssemblerNonlinear(p,e,t, area_vec, B, c, xi + k_n*k3)-epsAassembler(p,e,t,res, xi+k_n*k3, area_vec, B,c, beta,h))*(xi+k_n*k3));

        u0=xi;
        xi = u0 + k_n/6*(k1+2*k2+2*k3 + k4);
        xi(e(1,:)) = u_init(e(1,:));
        
        disp(i)

    end
    e = xi - u_init;
    L2E = sqrt(e'*M*e);

end

function [xi, u_init, L2E] = GFEM_RK4_solver(p,e,t, T, h_max, IC,field,k_n)
    [~,C,M, u_init, beta, ~, area_vec, B, c, ~, ~, ~] = assemble(p,e,t, IC,field);
    if (k_n == -1)
        CFL = 0.5;
        k_n = CFL*h_max/max(beta)*0.1;
    end
    k_n = T/ceil(T/k_n); 
    u0 = u_init;

    for i = 1:(T/k_n)

        k1 = M\(-convectionAssemblerNonlinear(p,e,t, area_vec,B, c, u0)*u0);
        k2 = M\(-convectionAssemblerNonlinear(p,e,t, area_vec,B, c, u0+k_n/2*k1)*(u0+k_n/2*k1));
        k3 = M\(-convectionAssemblerNonlinear(p,e,t, area_vec,B, c, u0+k_n/2*k2)*(u0+k_n/2*k2));
        k4 = M\(-convectionAssemblerNonlinear(p,e,t, area_vec,B, c, u0+k_n*k3)*(u0+k_n*k3));
        xi = u0 + k_n/6*(k1+2*k2+2*k3 + k4);
        xi(e(1,:)) = u_init(e(1,:));
        u0 = xi; 
        disp(i)


    end
    e = xi - u_init;
    L2E = sqrt(e'*M*e);


end



function [epsA, eps_vec] = epsAassembler(p,e,t, res, xi, area_vec, B, c, beta, h)
    N = size(p,2);
    epsA = sparse(N,N);
    eps_vec = zeros(N,1);
    Cvel = 0.25;
    Crv = 1.0;
    normalization = norm(xi-mean(xi), Inf);
    test_res = res/normalization;
    for K = 1:size(t,2)
       
        nodes = t(1:3, K);        
        area_K = area_vec(K);
        beta_K = beta(K);
 
        eps_vec(K) = min([Cvel*h(K)*beta_K, Crv*h(K)^2*norm(test_res(nodes),Inf)]);
        eps = eps_vec(K);
        
        epsAK = (B(:,K)*B(:,K)' + c(:,K)*c(:,K)')*area_K;
        epsA(nodes,nodes) = epsA(nodes,nodes) +eps*epsAK;
    end
    

end

function [C, delta_C, delta_S] = updateSUPGmatrices(p,e,t,delta, area_vec, B, c,u)
   
    N=size(p,2);
    N_t=size(t,2);
    delta_S=sparse(N,N);
    delta_C=sparse(N,N);
    C=sparse(N,N);
    for K=1:N_t
        nodes=t(1:3,K);
        [fx, fy] = convectionfieldNonlinear(u(nodes));
        cxmid = mean(fx);
        cymid = mean(fy);
        
        
        dSdK=delta(K)*(cxmid*B(:,K)+cymid*c(:,K))*(cxmid*B(:,K)+cymid*c(:,K))'*area_vec(K);
        dCK = delta(K)*ones(3,1)*(cxmid*B(:,K)+cymid*c(:,K))'*area_vec(K)/3;
        CK = ones(3,1)*(cxmid*B(:,K)+cymid*c(:,K))'*area_vec(K)/3;
        
        C(nodes,nodes)=C(nodes,nodes)+CK;
        delta_C(nodes,nodes)=delta_C(nodes,nodes)+dCK;
        delta_S(nodes,nodes)=delta_S(nodes,nodes)+dSdK;
    end
end


function C = convectionAssemblerNonlinear(p,e,t,area_vec, B,c,u)
    N=size(p,2);
    C=sparse(N,N);
    for K=1:size(t,2)
        nodes=t(1:3,K);
        [fx, fy] = convectionfieldNonlinear(u(nodes));
        cxmid = mean(fx);
        cymid = mean(fy);
        CK = ones(3,1)*(cxmid*B(:,K)+cymid*c(:,K))'*area_vec(K)/3;
        
        C(nodes,nodes)=C(nodes,nodes)+CK;
    end


end

function delta_C = deltaConvectionAssemblerNonlinear(p,e,t, delta, area_vec, B, c,u)
    N=size(p,2);
    
    delta_C=sparse(N,N);
    for K=1:size(t,2)
        nodes=t(1:3,K);         
        [fx, fy] = convectionfieldNonlinear(u(nodes));
        cxmid = mean(fx);
        cymid = mean(fy);  
        dCK = delta(K)*ones(3,1)*(cxmid*B(:,K)+cymid*c(:,K))'*area_vec(K)/3;
        
        delta_C(nodes,nodes)=delta_C(nodes,nodes)+dCK;
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

function delta_Sd = deltaSDAssemblerNonlinear(p,t, delta, area_vec, B, c, u)
    N=size(p,2);
    N_t=size(t,2);
    delta_Sd=sparse(N,N);
    for K=1:N_t
        nodes=t(1:3,K);
        [fx, fy] = convectionfieldNonlinear(u(nodes));
        cxmid = mean(fx);
        cymid = mean(fy);
        
        dSdK=delta(K)*(cxmid*B(:,K)+cymid*c(:,K))*(cxmid*B(:,K)+cymid*c(:,K))'*area_vec(K);
        delta_Sd(nodes,nodes)=delta_Sd(nodes,nodes)+dSdK;
    end
end


function [A,C,M, u_init, beta, h, area_vec, B, c, bx, by, delta] = assemble(p,e,t,fn, field)
    N=size(p,2);
    N_t = size(t,2);
    A=sparse(N,N);
    b = zeros(N,1);
    
    C=sparse(N,N);
    M=sparse(N,N);
    u_init = fn(p(1,:),p(2,:))';
    if (func2str(field) == "convectionfieldNonlinear")
        [bx, by] = field(u_init);
    else
        [bx, by] = field(p(1,:),p(2,:));
    end
    
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
        
        CK = ones(3,1)*(cxmid*B_K+cymid*c_K)'*area_K/3;
        MK = [2 1 1; 1 2 1; 1 1 2]/12*area_K;
       
        A(nodes,nodes) = A(nodes,nodes) +AK;
        
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

function out = shockwaveIC(x,y)
    r0 = 1;
    x0 = 0;
    y0 = 0;
    for i = 1:size(x,2)
        if (((x(i)-x0)^2+(y(i)-y0)^2) <= r0^2)
            out(i) = 14*pi/4;%ones(size(x,2));
        else
            out(i) = pi/4;%zeros(size(x,2));
        end
    end

end

function [bx, by] = convectionfield(x,y)
    bx = 2*pi*(-y);
    by = 2*pi*x;
end

function [bx, by] = convectionfieldNonlinear(u)
    bx = cos(u);
    by = -sin(u);
end


function [area,b,c] = HatGradients(x,y)
    area=polyarea(x,y);
    b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
    c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end


function r = Rectg(xmin,ymin,xmax,ymax)
    r=[2 xmin xmax ymin ymin 1 0;
    2 xmax xmax ymin ymax 1 0;
    2 xmax xmin ymax ymax 1 0;
    2 xmin xmin ymax ymin 1 0]';
end