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

U_0 = u_init(p(1,:),p(2,:));
U_prev = U_0;

%step once to calculate viscosity
A = (2/k_n)*M+C;
b = ((2/k_n)*M-C)*U_prev;
A(e(1,:),:) = I(e(1,:),:);
b(e(1,:)) = 0;
U = A\b;

t = k_n;
while t < T
    %Want to end at t = T
    %so will at final timestep take a smaller step
    %to reach T exactly
    if t+k_n > T
        k_n = T-t;
    end
    
    Rv = Rv_Assembler(p,tri,k_n,M,C,U,U_prev,bx,by,p_a_m);
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

function Rv = Rv_Assembler(p,t,k_n,M,C,U,U_prev,bx,by,pre_alloc_matrix)
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
