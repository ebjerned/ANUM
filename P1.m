clear;
close all;

g = @circleg;
h_max = 1/16;
T=1;
[p,e,t] = initmesh(g,'hmax',h_max);
[A,C,M,u_init,beta] = assemble(p,e,t);
CFL = 0.5;
k_n = CFL*h_max/beta;
k_n = T/ceil(T/k_n);

u0 = u_init;
I = eye(length(p));

[fx fy] = convectionfield(p(1,:), p(2,:));
quiver(p(1,:), p(2,:), fx,fy);
axis([-1.1 1.1 -1.1 1.1]);

figure;

pdeplot(p,e,t,'XYData',u_init, "ZData", u_init);
%caxis([-10 2]);
title("Linear Advection Equation at t = "+ 0+" , h_{max} = " + h_max);
xlabel("x");
ylabel("y");
%view(2);
for i=1:(T/k_n)
    b = (M-k_n/2*C)*u0;
    A = (M+k_n/2*C);
    
    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:))=0;
    xi = A\b; 
%     if (mod(i,floor(T/(10*k_n))) == 0 || i*k_n == 1)
      if (i*k_n == 1)
        figure;
        pdeplot(p,e,t,'XYData',xi, "ZData", xi);
        title("Linear Advection Equation at t = "+ i*k_n+" , h_{max} = " + h_max);
        xlabel("x");
        ylabel("y");
        %caxis([-10 2]);
        %view(2);         
     end
    u0 = xi;
end
e = xi - u_init;
L2E = sqrt(e'*M*e)
err_vec = [0.1340, 0.0530, 0.0117, 0.0060];
err_vec2 = [0.2143,0.1838, 0.1410, 0.1260];
figure;
h_steps = [1/4, 1/8, 1/16, 1/32];

coeff = polyfit(log(h_steps), log(err_vec),1);
coeff2 = polyfit(log(h_steps), log(err_vec2),1);
axes('XScale', 'log', 'YScale', 'log')
box on
hold;
grid;
loglog(h_steps, h_steps);
loglog(h_steps, err_vec);
loglog(h_steps, err_vec2);
loglog(h_steps, h_steps.^coeff(1).*exp(coeff(2)));
loglog(h_steps, h_steps.^coeff2(1).*exp(coeff2(2)));
legend(["\alpha = 1","1.1 Measured", "1.3 Measured", "1.1 \alpha = "+coeff(1), "1.3 \alpha= " + coeff2(1)]);
title("Convergence rate of different IC for GFEM");
xlabel("h_{max} [-]");
ylabel("L^2-norm error [-]");


function [A,C,M, u_init, beta] = assemble(p,e,t)
    N=size(p,2);
    A=sparse(N,N);
    b = zeros(N,1);
    C=sparse(N,N);
    M=sparse(N,N);
    u_init = f(p(1,:),p(2,:))';
    eps = zeros(N,1);
    Cvel = 0.25;
    Crv = 1.0;
    
    
    
    %Res = 1/k_n(u(end, :)-u(end-1,:)) + u_init*u(end,:);
    
    %Res = Res/;
    
    for K = 1:size(t,2)
        %page 243 in book
        nodes = t(1:3, K);
        x=p(1,nodes);
        y=p(2,nodes);
        [bx,by] = convectionfield(x,y);
        
        xc = mean(x);
        yc = mean(y);
        fc = f(xc, yc);
        cxmid = mean(bx);
        cymid = mean(by);
        
        beta_K = max(sqrt(bx.^2 + by.^2));
        beta(K)= beta_K;
        %Res_K = ;
        %eps(K) = min(Cvel*h*beta_K, Crv*h^2*Res_K);
        
        
        [area_K, B,c] = HatGradients(x,y);
        AK = (B*B' + c*c')*area_K; % Need to multiply with residual viscosity
        bK = fc*area_K/3;
        
        CK = ones(3,1)*(cxmid*B+cymid*c)'*area_K/3;
        
        MK = [2 1 1; 1 2 1; 1 1 2]/12*area_K;
        A(nodes,nodes) = A(nodes,nodes) +AK;
        b(nodes)= b(nodes) + bK;
        C(nodes,nodes) = C(nodes,nodes) +CK;
        M(nodes,nodes) = M(nodes,nodes) +MK;
    
        
    
    end
    beta = max(beta_K);
end

%function eps = res_viscosity(

function out = f(x,y)
    r0 = 0.25;
    x0 = 0.3;
    y0 = 0;
    out = 0.5*(1-tanh((((x-x0).^2+(y-y0).^2))./r0^2-1));
%     for i = 1:size(x,2)
%         if (((x(i)-x0)^2+(y(i)-y0)^2) <= r0^2)
%            out(i) = 1;%ones(size(x,2));
%         else
%             out(i) = 0;%zeros(size(x,2));
%         end
%     end
    
    
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