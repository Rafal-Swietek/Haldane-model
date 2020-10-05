function Berry_kpath_kagome(k_path,grid,parameters,a1,a2)
%Calculating berry phase
t1 = parameters(1); L1 = parameters(2); 
t2 = parameters(3); L2 = parameters(4);
F_1 = zeros(1,length(k_path));
F_2 = zeros(1,length(k_path));
F_3 = zeros(1,length(k_path));
Chern1 = 0; Chern2 = 0;
dk = 2*pi/3/grid;
for ii=1:length(k_path)
    %first point (bottom left)
        k_x = k_path(1,ii)-dk/2;
        k_y = k_path(2,ii)-dk/2;
            k1 = k_x*a1(1)+k_y*a1(2);    k2 = k_x*a2(1)+k_y*a2(2);
            u = exp(i*k1);
            v = exp(i*k2); %conj(u) = 1/u, because its exponential
            
                %Nearest Neighbours hamiltonian
                H_NN(1,2) = -(t1+i*L1)*(1+u); H_NN(2,1) = conj(H_NN(1,2));
                H_NN(1,3) = -(t1-i*L1)*(1+v); H_NN(3,1) = conj(H_NN(1,3));
                H_NN(2,3) = -(t1+i*L1)*(1+v/u); H_NN(3,2) = conj(H_NN(2,3));
                %Next Nearest Neighbours hamiltonian
                H_NNN(1,2) = -(t2-i*L2)*(v+u/v); H_NNN(2,1) = conj(H_NNN(1,2));
                H_NNN(1,3) = -(t2+i*L2)*(u+v/u); H_NNN(3,1) = conj(H_NNN(1,3));
                H_NNN(2,3) = -(t2-i*L2)*(v+1/u); H_NNN(3,2) = conj(H_NNN(2,3));
            %Hamiltonian + eigenvalue
            H = H_NN + H_NNN; 
        [Psi1,~] = eig(H);
                Psi1(1,:) = Psi1(1,:)*exp(i*(k_x*0+k_y*0)); %bloch function stay in matrix Psi
                Psi1(2,:) = Psi1(2,:)*exp(i*k1/2);
                Psi1(3,:) = Psi1(3,:)*exp(i*k2/2);
          
    %second point (bottom right)
        k_x = k_path(1,ii)+dk/2;
        k_y = k_path(2,ii)-dk/2;
            k1 = k_x*a1(1)+k_y*a1(2);    k2 = k_x*a2(1)+k_y*a2(2);
            u = exp(i*k1);
            v = exp(i*k2); %conj(u) = 1/u, because its exponential
            
                %Nearest Neighbours hamiltonian
                H_NN(1,2) = -(t1+i*L1)*(1+u); H_NN(2,1) = conj(H_NN(1,2));
                H_NN(1,3) = -(t1-i*L1)*(1+v); H_NN(3,1) = conj(H_NN(1,3));
                H_NN(2,3) = -(t1+i*L1)*(1+v/u); H_NN(3,2) = conj(H_NN(2,3));
                %Next Nearest Neighbours hamiltonian
                H_NNN(1,2) = -(t2-i*L2)*(v+u/v); H_NNN(2,1) = conj(H_NNN(1,2));
                H_NNN(1,3) = -(t2+i*L2)*(u+v/u); H_NNN(3,1) = conj(H_NNN(1,3));
                H_NNN(2,3) = -(t2-i*L2)*(v+1/u); H_NNN(3,2) = conj(H_NNN(2,3));
            %Hamiltonian + eigenvalue
            H = H_NN + H_NNN;
        [Psi2,~] = eig(H);
                Psi2(1,:) = Psi2(1,:)*exp(i*(k_x*0+k_y*0)); %bloch function stay in matrix Psi
                Psi2(2,:) = Psi2(2,:)*exp(i*k1/2);
                Psi2(3,:) = Psi2(3,:)*exp(i*k2/2);
          
    %third point (top right)
        k_x = k_path(1,ii)+dk/2;
        k_y = k_path(2,ii)+dk/2;
            k1 = k_x*a1(1)+k_y*a1(2);    k2 = k_x*a2(1)+k_y*a2(2);
            u = exp(i*k1);
            v = exp(i*k2); %conj(u) = 1/u, because its exponential
            
                %Nearest Neighbours hamiltonian
                H_NN(1,2) = -(t1+i*L1)*(1+u); H_NN(2,1) = conj(H_NN(1,2));
                H_NN(1,3) = -(t1-i*L1)*(1+v); H_NN(3,1) = conj(H_NN(1,3));
                H_NN(2,3) = -(t1+i*L1)*(1+v/u); H_NN(3,2) = conj(H_NN(2,3));
                %Next Nearest Neighbours hamiltonian
                H_NNN(1,2) = -(t2-i*L2)*(v+u/v); H_NNN(2,1) = conj(H_NNN(1,2));
                H_NNN(1,3) = -(t2+i*L2)*(u+v/u); H_NNN(3,1) = conj(H_NNN(1,3));
                H_NNN(2,3) = -(t2-i*L2)*(v+1/u); H_NNN(3,2) = conj(H_NNN(2,3));
            %Hamiltonian + eigenvalue
            H = H_NN + H_NNN;
        [Psi3,~] = eig(H);
                Psi3(1,:) = Psi3(1,:)*exp(i*(k_x*0+k_y*0)); %bloch function stay in matrix Psi
                Psi3(2,:) = Psi3(2,:)*exp(i*k1/2);
                Psi3(3,:) = Psi3(3,:)*exp(i*k2/2);
          
    %fourth point (top left)
        k_x = k_path(1,ii)-dk/2;
        k_y = k_path(2,ii)+dk/2;
            k1 = k_x*a1(1)+k_y*a1(2);    k2 = k_x*a2(1)+k_y*a2(2);
            u = exp(i*k1);
            v = exp(i*k2); %conj(u) = 1/u, because its exponential
            
                %Nearest Neighbours hamiltonian
                H_NN(1,2) = -(t1+i*L1)*(1+u); H_NN(2,1) = conj(H_NN(1,2));
                H_NN(1,3) = -(t1-i*L1)*(1+v); H_NN(3,1) = conj(H_NN(1,3));
                H_NN(2,3) = -(t1+i*L1)*(1+v/u); H_NN(3,2) = conj(H_NN(2,3));
                %Next Nearest Neighbours hamiltonian
                H_NNN(1,2) = -(t2-i*L2)*(v+u/v); H_NNN(2,1) = conj(H_NNN(1,2));
                H_NNN(1,3) = -(t2+i*L2)*(u+v/u); H_NNN(3,1) = conj(H_NNN(1,3));
                H_NNN(2,3) = -(t2-i*L2)*(v+1/u); H_NNN(3,2) = conj(H_NNN(2,3));
            %Hamiltonian + eigenvalue
            H = H_NN + H_NNN;
        [Psi4,~] = eig(H);
                Psi4(1,:) = Psi4(1,:)*exp(i*(k_x*0+k_y*0)); %bloch function stay in matrix Psi
                Psi4(2,:) = Psi4(2,:)*exp(i*k1/2);
                Psi4(3,:) = Psi4(3,:)*exp(i*k2/2);
     
          
        x12 = dot(Psi1(:,1),Psi2(:,1));
        x23 = dot(Psi2(:,1),Psi3(:,1));
        x34 = dot(Psi3(:,1),Psi4(:,1));
        x41 = dot(Psi4(:,1),Psi1(:,1));
        berry_1 = -angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
        F_1(ii) = berry_1/(dk^2);
        
        x12 = dot(Psi1(:,2),Psi2(:,2));
        x23 = dot(Psi2(:,2),Psi3(:,2));
        x34 = dot(Psi3(:,2),Psi4(:,2));
        x41 = dot(Psi4(:,2),Psi1(:,2));
        berry_2 = -angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
        F_2(ii) = berry_2/(dk^2);
        
        x12 = dot(Psi1(:,3),Psi2(:,3));
        x23 = dot(Psi2(:,3),Psi3(:,3));
        x34 = dot(Psi3(:,3),Psi4(:,3));
        x41 = dot(Psi4(:,3),Psi1(:,3));
        berry_3 = -angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
        F_3(ii) = berry_3/(dk^2);
        % minus sign by angle, because other sided loop around pint 
    
end
%%
% %Berry curvature plot

figure(16);
plot(1:length(k_path(1,:)),F_1,'b-',1:length(k_path(1,:)),F_2,'r-',1:length(k_path(1,:)),F_3,'k-');
title(sprintf('Berry curvature along k-path: \\Gamma-K-M-\\Gamma \n using parameters: t1 =%1.1f, t2 = %1.2f, \\lambda1 = %1.2f and \\lambda2 = %1.2f ',t1,t2,L1,L2));
xticks([0 length(k_path)/3 1.5*length(k_path)/3 2*length(k_path)/3 length(k_path)]);   
xticklabels({'\Gamma','K','M','K"','\Gamma'});
axis([0 length(k_path) 1.01*min([min(F_1),min(F_2),min(F_3)]) 1.01*max([max(F_1),max(F_2),max(F_3)])]);
ylabel(sprintf('\\Omega_z')); 
legend('3rd band','2nd band', '1st band');

end

