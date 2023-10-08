[u,p,e,tri] = PDE_Solve_SUPG(1/4,1,1e-3,1e-1,@shockwaveIC);

function [u,p,e,tri] = PDE_Solve_SUPG(hmax,T,k_n,pic_TOL,u_init)
g = Rectg(-2,2,-2.5,1.5);
[p,e,tri] = initmesh(g,'hmax',hmax);
u_0 = u_init(p(1,:),p(2,:))'; 
u_old = u_0;
u = u_old;

%assemble matrices
[bx,by] = f_prim(u);
p_a_m = pre_alloc_matrix(p,tri); %[area,b,c,h_K]
delta_pam = delta_add_p_a_m(tri,bx,by,p_a_m);
M = Mass_assembler_2D(p,tri,delta_pam);

t = 0;
while t < T
    if t+k_n > T
        k_n = T-t;
    end
    disp("-----------------")
    disp("Solver: SUPG")
    disp("Time: " + num2str(t+k_n))
    u_old_tmp = u;%to save u_old to next step
    pic_err = pic_TOL+1;
    while pic_err > pic_TOL
        u_tmp = u;
        [bx,by] = f_prim(u_tmp);
        delta_pam = delta_add_p_a_m(tri,bx,by,p_a_m);
        C = ConvectionAssembler2D(p,tri,bx,by,delta_pam);
        H = H_assembler(p,tri,bx,by,delta_pam);
        Sd = SD_assembler(p,tri,bx,by,delta_pam);
        
        A = (2/k_n)*(M+H)+C+Sd;
        b = ((2/k_n)*(M+H)-C-Sd)*u_old;
        u = A\b;
        
        %apply BC
        u(e(1,:)) = u_0(e(1,:));
        
        pic_err = norm(u_tmp-u);
        disp("Picard error: " + num2str(pic_err))
    end
    u_old = u_old_tmp;
    t = t+k_n;
end

function r = Rectg(xmin,xmax,ymin,ymax)
r=[2 xmin xmax ymin ymin 1 0;
2 xmax xmax ymin ymax 1 0;
2 xmax xmin ymax ymax 1 0;
2 xmin xmin ymax ymin 1 0]';
end

function [bx,by] = f_prim(u)
    bx = cos(u);
    by = -sin(u);
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

function out = delta_add_p_a_m(t,bx,by,pre_alloc_matrix)
   nt = size(t,2);
   out = pre_alloc_matrix;
   for k = 1:nt
       h_k = out(k,8);
       
       %calculated delta
       loc2glb = t(1:3,k);
       delta_K = 0.5*h_k/norm([bx(loc2glb),by(loc2glb)]);
       
       out(k,8) = delta_K;
   end
end

function M = Mass_assembler_2D(p,t,pre_alloc_matrix)
%MASS_ASSEMBLER_2D
%input is point matrix p and connectivity matrix t given by
%initmesh
%also matrix consisting of hat functions
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

function Sd = SD_assembler(p,t,bx,by,pre_alloc_matrix)
np=size(p,2);
nt=size(t,2);
Sd=sparse(np,np);
for k=1:nt
    loc2glb=t(1:3,k);
    area = pre_alloc_matrix(k,1);
    b_SD = pre_alloc_matrix(k,2:4)';
    c = pre_alloc_matrix(k,5:7)';
    delta_K = pre_alloc_matrix(k,8);
    bxmid=mean(bx(loc2glb));
    bymid=mean(by(loc2glb));
    SdK=(bxmid*b_SD+bymid*c)*(bxmid*b_SD+bymid*c)'*area;
    
    Sd(loc2glb,loc2glb)=Sd(loc2glb,loc2glb)+delta_K*SdK;
end
end

function H = H_assembler(p,t,bx,by,pre_alloc_matrix)
np=size(p,2);
nt=size(t,2);
H=sparse(np,np);
for k=1:nt
    loc2glb=t(1:3,k);
    area = pre_alloc_matrix(k,1);
    b_H = pre_alloc_matrix(k,2:4)';
    c = pre_alloc_matrix(k,5:7)';
    delta_K = pre_alloc_matrix(k,8);
    bxmid=mean(bx(loc2glb));
    bymid=mean(by(loc2glb));
    SdK=ones(3,1)*(bxmid*b_H+bymid*c)'*area/3;
    
    H(loc2glb,loc2glb)=H(loc2glb,loc2glb)+delta_K*SdK;
end
H = H';%H = C'
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
