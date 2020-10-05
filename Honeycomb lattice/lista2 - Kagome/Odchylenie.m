function Odchylenie(parameters,kx,ky,a1,a2)
t1 = parameters(1); L1 = parameters(2); 
t2 = parameters(3); L2 = parameters(4);
    KX = 4*pi/9*sqrt(3);
    Kax = [KX,KX/2,-KX/2,-KX,-KX/2,KX/2,KX]/norm(a1);
    Kay = [0,2*pi/3,2*pi/3,0,-2*pi/3,-2*pi/3,0]/norm(a2);
    F_1 = zeros(length(kx),length(ky));
    F_2 = zeros(length(kx),length(ky));
    F_3 = zeros(length(kx),length(ky));
    Chern1 = 0; Chern2 = 0; Chern3 = 0;
    step = 20;
    %L1 = 0.5;
    lambda1 = -0.2:1/step:0.2;
    lambda2 = -0.2:1/step:0.2;
    Odchylenie1 = zeros(length(lambda1),length(lambda2));
    Odchylenie2 = zeros(length(lambda1),length(lambda2));
    Odchylenie3 = zeros(length(lambda1),length(lambda2));
for ID2 = 1:length(lambda1)
    L1 = lambda1(ID2)
 for ID = 1:length(lambda1)
    L2 = lambda2(ID);
    Fav1 = 0; Fav2 = 0; Fav3 = 0;
    Favsq1 = 0; Favsq2 = 0; Favsq3 = 0;
    iter = 0;   
    for jj=1:length(ky)-1
       for ii=1:length(kx)-1
        if((abs(kx(ii))<4*pi/3&&(abs(ky(jj))<2*pi/3*sqrt(3))))
         if((ky(jj)<(-sqrt(3)*kx(ii)+4*pi/sqrt(3)))&&(ky(jj)>(sqrt(3)*kx(ii)-4*pi)))
          if((ky(jj)<(sqrt(3)*kx(ii)+4*pi/sqrt(3)))&&(ky(jj)>(-sqrt(3)*kx(ii)-4*pi/sqrt(3))))
                % first point
                    k_x = kx(ii);
                    k_y = ky(jj);
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
                        Psi1(2,:) = Psi1(2,:)*exp(i*k1);
                        Psi1(3,:) = Psi1(3,:)*exp(i*k2);
                 %second point
                    k_x = kx(ii+1);
                    k_y = ky(jj);
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
                        Psi2(2,:) = Psi2(2,:)*exp(i*k1);
                        Psi2(3,:) = Psi2(3,:)*exp(i*k2);

                  %third point
                    k_x = kx(ii+1);
                    k_y = ky(jj+1);
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
                        Psi3(2,:) = Psi3(2,:)*exp(i*k1);
                        Psi3(3,:) = Psi3(3,:)*exp(i*k2);

                    %fourth point
                    k_x = kx(ii);
                    k_y = ky(jj+1);
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
                        Psi4(2,:) = Psi4(2,:)*exp(i*k1);
                        Psi4(3,:) = Psi4(3,:)*exp(i*k2);

                            x =    dot(Psi1(:,1),Psi2(:,1))/abs(dot(Psi1(:,1),Psi2(:,1)));
                            x = x *dot(Psi2(:,1),Psi3(:,1))/abs(dot(Psi2(:,1),Psi3(:,1)));
                            x = x *dot(Psi3(:,1),Psi4(:,1))/abs(dot(Psi3(:,1),Psi4(:,1)));
                            x = x *dot(Psi4(:,1),Psi1(:,1))/abs(dot(Psi4(:,1),Psi1(:,1)));
                            berry_1 = angle(x);

                            x =    dot(Psi1(:,2),Psi2(:,2))/abs(dot(Psi1(:,2),Psi2(:,2)));
                            x = x *dot(Psi2(:,2),Psi3(:,2))/abs(dot(Psi2(:,2),Psi3(:,2)));
                            x = x *dot(Psi3(:,2),Psi4(:,2))/abs(dot(Psi3(:,2),Psi4(:,2)));
                            x = x *dot(Psi4(:,2),Psi1(:,2))/abs(dot(Psi4(:,2),Psi1(:,2)));
                            berry_2 = angle(x);

                            x =    dot(Psi1(:,3),Psi2(:,3))/abs(dot(Psi1(:,3),Psi2(:,3)));
                            x = x *dot(Psi2(:,3),Psi3(:,3))/abs(dot(Psi2(:,3),Psi3(:,3)));
                            x = x *dot(Psi3(:,3),Psi4(:,3))/abs(dot(Psi3(:,3),Psi4(:,3)));
                            x = x *dot(Psi4(:,3),Psi1(:,3))/abs(dot(Psi4(:,3),Psi1(:,3)));
                            berry_3 = angle(x);

                            d2k = ( (kx(ii)-kx(ii+1))^2 + (ky(jj)-ky(jj))^2 )^(1/2)*( (kx(ii+1)-kx(ii))^2 + (ky(jj+1)-ky(jj+1))^2 )^(1/2);

                            F_1(ii,jj) = berry_1/d2k; 
                            F_2(ii,jj) = berry_2/d2k; 
                            F_3(ii,jj) = berry_3/d2k;
                            
                              Fav1 = Fav1 + F_1(ii,jj);
                              Fav2 = Fav2 + F_2(ii,jj);
                              Fav3 = Fav3 + F_3(ii,jj);
                              Favsq1 = Favsq1 + (F_1(ii,jj))^2;
                              Favsq2 = Favsq2 + (F_2(ii,jj))^2;
                              Favsq3 = Favsq3 + (F_3(ii,jj))^2;
                              iter = iter+1;
            end
           end    
          end
        end
    end
    Odchylenie1(ID2,ID) = sqrt(-(Fav1/iter)^2 + Favsq1/iter);
    Odchylenie2(ID2,ID) = sqrt(-(Fav2/iter)^2 + Favsq2/iter);
    Odchylenie3(ID2,ID) = sqrt(-(Fav3/iter)^2 + Favsq3/iter);
 end
end
% plot(lambda1,Odchylenie1,'k-',lambda1,Odchylenie2,'b-',lambda1,Odchylenie3,'r-');
figure(131);
mesh(lambda1,lambda2,Odchylenie1);
title(sprintf('Standard deviaton of Berry curvature for 1st band in function of parameters \\lambda_2 & \\lambda_1\n for other parameters: t_1 = %0.2f and t_2 = %0.2f',t1,t2));
xlabel(sprintf('\\lambda_1'));
ylabel(sprintf('\\lambda_2'));
zlabel(sprintf('\\sigma(\\Omega_z)'));

figure(4234);
mesh(lambda1,lambda2,Odchylenie2);
title(sprintf('Standard deviaton of Berry curvature for 2nd band in function of parameters \\lambda_2 & \\lambda_1\n for other parameters: t_1 = %0.2f and t_2 = %0.2f',t1,t2));
xlabel(sprintf('\\lambda_1'));
ylabel(sprintf('\\lambda_2'));
zlabel(sprintf('\\sigma(\\Omega_z)'));

figure(185);
mesh(lambda1,lambda2,Odchylenie3);
title(sprintf('Standard deviaton of Berry curvature for 3rd band in function of parameters \\lambda_2 & \\lambda_1\n for other parameters: t_1 = %0.2f and t_2 = %0.2f',t1,t2));
xlabel(sprintf('\\lambda_1'));
ylabel(sprintf('\\lambda_2'));
zlabel(sprintf('\\sigma(\\Omega_z)'));
%legend('3rd band','2nd band', '1st band');

    
end



