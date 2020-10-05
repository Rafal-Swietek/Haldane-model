function Berry_k_path(k_path,parameters,grid,a1,a2)
%Calculating berry phase
t = parameters(1); V = parameters(2); lambda = parameters(3);
F_C = zeros(1,length(k_path));
F_V = zeros(1,length(k_path));
Chern1 = 0; Chern2 = 0;
dk = 4*pi*sqrt(3)/9/grid;
for ii=1:length(k_path)
    %first point (bottom left)
        k_x = k_path(1,ii)-dk/2;
        k_y = k_path(2,ii)-dk/2;
        k1 = k_x*a1(1)+k_y*a1(2);
        k2 = k_x*a2(1)+k_y*a2(2);
        k12 = k_x*(a1(1)-a2(1))+k_y*(a1(2)-a2(2));
        %D-vector
        d_x = t*(1 + cos(k1) + cos(k2));
        d_y = t*(sin(k1) + sin(k2));
        d_z = 2*V - 4*lambda*( sin(k1) - sin(k2) - sin(k12));
        %Hamiltonian + eigenvalue
        H = [d_z, d_x-i*d_y; d_x+i*d_y, -d_z];
        [Psi1,~] = eig(H);
          
    %second point (bottom right)
        k_x = k_path(1,ii)+dk/2;
        k_y = k_path(2,ii)-dk/2;
        k1 = k_x*a1(1)+k_y*a1(2);
        k2 = k_x*a2(1)+k_y*a2(2);
        k12 = k_x*(a1(1)-a2(1))+k_y*(a1(2)-a2(2));
        %D-vector
        d_x = t*(1 + cos(k1) + cos(k2));
        d_y = t*(sin(k1) + sin(k2));
        d_z = 2*V - 4*lambda*( sin(k1) - sin(k2) - sin(k12));
        %Hamiltonian + eigenvalue
        H = [d_z, d_x-i*d_y; d_x+i*d_y, -d_z];
        [Psi2,~] = eig(H);
          
    %third point (top right)
        k_x = k_path(1,ii)+dk/2;
        k_y = k_path(2,ii)+dk/2;
        k1 = k_x*a1(1)+k_y*a1(2);
        k2 = k_x*a2(1)+k_y*a2(2);
        k12 = k_x*(a1(1)-a2(1))+k_y*(a1(2)-a2(2));
        %D-vector
        d_x = t*(1 + cos(k1) + cos(k2));
        d_y = t*(sin(k1) + sin(k2));
        d_z = 2*V - 4*lambda*( sin(k1) - sin(k2) - sin(k12));
        %Hamiltonian + eigenvalue
        H = [d_z, d_x-i*d_y; d_x+i*d_y, -d_z];
        [Psi3,~] = eig(H);
          
    %fourth point (top left)
        k_x = k_path(1,ii)-dk/2;
        k_y = k_path(2,ii)+dk/2;
        k1 = k_x*a1(1)+k_y*a1(2);
        k2 = k_x*a2(1)+k_y*a2(2);
        k12 = k_x*(a1(1)-a2(1))+k_y*(a1(2)-a2(2));
        %D-vector
        d_x = t*(1 + cos(k1) + cos(k2));
        d_y = t*(sin(k1) + sin(k2));
        d_z = 2*V - 4*lambda*( sin(k1) - sin(k2) - sin(k12));
        %Hamiltonian + eigenvalue
        H = [d_z, d_x-i*d_y; d_x+i*d_y, -d_z];
        [Psi4,~] = eig(H);
     
          
        x12 = dot(Psi1(:,1),Psi2(:,1));
        x23 = dot(Psi2(:,1),Psi3(:,1));
        x34 = dot(Psi3(:,1),Psi4(:,1));
        x41 = dot(Psi4(:,1),Psi1(:,1));
        berry_C = angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
        F_C(ii) = berry_C/(dk^2);
        
        x12 = dot(Psi1(:,2),Psi2(:,2));
        x23 = dot(Psi2(:,2),Psi3(:,2));
        x34 = dot(Psi3(:,2),Psi4(:,2));
        x41 = dot(Psi4(:,2),Psi1(:,2));
        berry_V = angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
        F_V(ii) = berry_V/(dk^2);
    
end
%%
%Berry phase&curvature plot
chern_cond = Chern1/(2*pi);
chern_val = Chern2/(2*pi);

figure(14);
plot(1:length(k_path(1,:)),F_C,'k-');
title(sprintf('Berry curvature for conducting band \n using parameters: V =%1.2f, t = %1.1f and \\lambda = %1.2f ',V,t,lambda));
xticks([0 length(k_path)/3 1.5*length(k_path)/3 2*length(k_path)/3 length(k_path)]);   
xticklabels({'\Gamma','K','M','K"','\Gamma'});
axis([0 length(k_path) min(F_C) max(F_C)]);
ylabel(sprintf('\\Omega_z'));    

figure(15);
plot(1:length(k_path(1,:)),F_V,'k-');
title(sprintf('Berry curvature for valence band \n using parameters: V =%1.2f, t = %1.1f and \\lambda = %1.2f ',V,t,lambda));
xticks([0 length(k_path)/3 1.5*length(k_path)/3 2*length(k_path)/3 length(k_path)]);   
xticklabels({'\Gamma','K','M','K"','\Gamma'});
axis([0 length(k_path) min(F_V) max(F_V)]);
ylabel(sprintf('\\Omega_z')); 

end

