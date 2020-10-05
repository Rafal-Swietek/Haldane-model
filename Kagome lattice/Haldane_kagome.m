%K-space base vectors
a = 1; %lattice constant

a1 = [a, 0];    %Bravais lattice vectors
a2 = [a/2, a*sqrt(3)/2];

syms b_1x b_1y b_2x b_2y;
equations1 = [b_1x*a1(1) + b_1y*a1(2)==2*pi, b_1x*a2(1) + b_1y*a2(2)==0   ];
equations2 = [b_2x*a1(1) + b_2y*a1(2)== 0  , b_2x*a2(1) + b_2y*a2(2)==2*pi];
sol =   solve(equations1, [b_1x b_1y]); %solving for b1
solut = solve(equations2, [b_2x b_2y]); %solving for b2
b1 = [sol.b_1x sol.b_1y]; %assigning solution to vector
b2 = [solut.b_2x solut.b_2y];  %-||-
%%
% Finding K & M points
KK = [b1; b2; b1+b2; -b1; -b2; -b1-b2]; %Unit cell
r = norm(b1)/3^0.5; %radius of circle on Wigner-Seitz cell
alfa = pi/3;  
Rot = [cos(alfa) -sin(alfa); sin(alfa) cos(alfa)]; %matrix rotating around 60 degrees (pi/3)
K1 = [r 0]';
K = [ K1'; (Rot*K1)'; (Rot*(Rot*K1))'; -K1'; -(Rot*K1)'; -(Rot*(Rot*K1))' ]; %high symmetry K-points
M1 = (b1+b2)'/2;
M = [M1'; (Rot*M1)'; (Rot*(Rot*M1))'; -M1'; -(Rot*M1)'; -(Rot*(Rot*M1))']; %all M points
%%
%Energy eigenvalue calculation

                        t1 = 1; %hopping integral
                        L1 = 0.6; %complex hopping for NN
                        t2 = -0.3; % NNN hopping
                        L2 = 0.0; % next nearest neighbours hopping
                        
    grid_k = 20;
    %Array initializing
    kx = -pi/norm(a1):pi/norm(a1)/grid_k:pi/norm(a1); %kx grid
    ky = -pi/norm(a2):pi/norm(a2)/grid_k:pi/norm(a2); %ky grid
        E1 = zeros(length(kx),length(ky));
        E2 = zeros(length(kx),length(ky));
        E3 = zeros(length(kx),length(ky));
        E_kx = zeros(3,length(kx)); %energy along ky=0
        E_ky = zeros(3,length(ky)); %energy along kx=0
    %Simulation
   U = zeros(3,3,length(kx),length(ky)); %tensor for eigenvectors
   H = zeros(3,3);
   H_NN = H; H_NNN = H;
    for jj = 1:length(kx)
        k_y = ky(jj);
        for ii = 1:length(ky)
            k_x = kx(ii);
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
            [Psi,energy] = eig(H);
                Psi(1,:) = Psi(1,:)*exp(i*(k_x*0+k_y*0)); %bloch function stay in matrix Psi
                Psi(2,:) = Psi(2,:)*exp(i*k1/2);
                Psi(3,:) = Psi(3,:)*exp(i*k2/2);
            U(:,:,ii,jj) = Psi; % 4-tensor
            E1(ii,jj) = energy(1,1);
            E2(ii,jj) = energy(2,2); 
            E3(ii,jj) = energy(3,3);
%             if(k_y==0)
%                 E_kx(1,ii) = E1(ii,jj);
%                 E_kx(2,ii) = E2(ii,jj);
%                 E_kx(3,ii) = E3(ii,jj);
%             end
%             if(k_x==0)
%                 E_ky(1,jj) = E1(ii,jj);
%                 E_ky(2,jj) = E2(ii,jj);
%                 E_ky(3,jj) = E3(ii,jj);
%             end
        end
    end
%%
%Defining k-path (œcie¿ka k) & Energy along k-path
grid = 100; %k_path grid
%k_path = K_gridGMKG(b1,b2,a1,a2,K,M,grid); %returns a vector with pints on kpath
% Chern = zeros(3,2*step+1);
% F_average = zeros(3,2*step+1);
parameters = [t1, L1, t2, L2];
%Energy_k_path_kagome(k_path,parameters,a1,a2);
%Berry_kagome(U,kx,ky,parameters,a1,a2);
%Berry_kpath_kagome(k_path,grid,parameters,a1,a2);
Odchylenie(parameters,kx,ky,a1,a2);

%     figure(111);
%     plot(1:length(Chern(1,:)),Chern(1,:),1:length(Chern(1,:)),Chern(2,:),1:length(Chern(1,:)),Chern(3,:));
%     title(sprintf('Chern number vs \\lambda_1'));
%     figure(123);
%     plot(1:length(F_average(1,:)),F_average(1,:),1:length(F_average(1,:)),F_average(2,:),1:length(F_average(1,:)),F_average(3,:));
%     title(sprintf('Averaged berry curvature vs \\lambda_1'));
%%
% 
%Plotting 3D plot E(kx,ky) and 2D plots (E(kx,0) & E(0,ky)
% figure(2); % 3D plot
% Colour = gradient(E2-E1);
% meshc(kx,ky,E1); hold on
% meshc(kx,ky,E2); hold on
% meshc(kx,ky,E3);
% xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
% xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% yticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
% yticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% xlabel(sprintf('k_x [\\pi/a]')); ylabel(sprintf('k_y [\\pi/a]')); zlabel('E [eV]');
% title(sprintf('Energy spectrum \n using parameters: t1 =%1.1f, t2 = %1.2f, \\lambda1 = %1.2f and \\lambda2 = %1.2f ',t1,t2,L1,L2))

% figure(3); %along x axis -> ky = 0
% plot(kx,E_kx(1,:),'k-', kx, E_kx(2,:),'k-',kx, E_kx(3,:),'k-');
% xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
% xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% xlabel(sprintf('k_x [\\pi/a]')); ylabel('E [eV]');
% title(sprintf('Energy spectrum along k_y = 0 \n using parameters: t1 =%1.1f, t2 = %1.2f, \\lambda1 = %1.2f and \\lambda2 = %1.2f ',t1,t2,L1,L2))
% 
% figure(4); % along y-axis  -> kx = 0
% plot(ky,E_ky(1,:),'k-', ky, E_ky(2,:),'k-',ky, E_ky(3,:),'k-');
% xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
% xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% xlabel(sprintf('k_y [\\pi/a]')); ylabel('E [eV]');
% title(sprintf('Energy spectrum along k_x = 0 \n using parameters: t1 =%1.1f, t2 = %1.2f, \\lambda1 = %1.2f and \\lambda2 = %1.2f ',t1,t2,L1,L2))
% 



