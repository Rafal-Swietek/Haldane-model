%K-space base vectors
a1 = 0.5*[-3^(1/2), 3]; %Bravais lattice vectors
a2 = 0.5*[3^(1/2), 3];
syms b_1x b_1y b_2x b_2y;
equations1 = [b_1x*a1(1) + b_1y*a1(2)==2*pi, b_1x*a2(1) + b_1y*a2(2)==0   ];
equations2 = [b_2x*a1(1) + b_2y*a1(2)== 0  , b_2x*a2(1) + b_2y*a2(2)==2*pi];
sol =   solve(equations1, [b_1x b_1y]); %solving for b1
solut = solve(equations2, [b_2x b_2y]); %solving for b2
b1 = [sol.b_1x sol.b_1y]; %assigning solution to vector
b2 = [solut.b_2x solut.b_2y];%-||-
%%
% Finding K & M points
KK = [b1; b2; b1+b2; -b1; -b2; -b1-b2]; %Unit cell
r = norm(b1)/3^0.5; %radius of circle on Wigner-Seitz cell
alfa = pi/3;  
Rot = [cos(alfa) -sin(alfa); sin(alfa) cos(alfa)]; %matrix rotating around 60 degrees (pi/3)
K1 = [r 0]';
K = [ K1'; (Rot*K1)'; (Rot*(Rot*K1))'; -K1'; -(Rot*K1)'; -(Rot*(Rot*K1))' ]; %high symmetry K-points
M1 = b2'/2;
M = [M1'; (Rot*M1)'; (Rot*(Rot*M1))'; -M1'; -(Rot*M1)'; -(Rot*(Rot*M1))']; %all M points
%%
%Energy eigenvalue calculation

                        t = 1; %hopping integral
                        V = 0.1; %staggered potential
                        lambda = -0.05; % next nearest neighbours hopping
                        
    grid_k = 50;
    %Array initializing
    kx = -2*pi/norm(a1):pi/norm(a1)/grid_k:2*pi/norm(a1); %kx grid
    ky = -2*pi/norm(a2):pi/norm(a2)/grid_k:2*pi/norm(a2); %ky grid
        E1 = zeros(length(kx),length(ky));
        E2 = zeros(length(kx),length(ky));
        d_x = zeros(length(kx),length(ky));
        d_y = zeros(length(kx),length(ky));
        d_z = zeros(length(kx),length(ky));
        E_kx = zeros(2,length(kx)); %energy along ky=0
        E_ky = zeros(2,length(ky)); %energy along kx=0
    %Simulation
    U = zeros(2,2,length(kx),length(ky)); %tensor for eigenvectors
    phase = zeros(length(kx),length(ky));
    for jj = 1:length(ky)
        k_y = ky(jj);
        for ii = 1:length(kx)
            k_x = kx(ii);
            k1 = k_x*a1(1)+k_y*a1(2);
            k2 = k_x*a2(1)+k_y*a2(2);
            k12 = k_x*(a1(1)-a2(1))+k_y*(a1(2)-a2(2));
            %D-vector
            d_x(ii,jj) = t*(1 + cos(k1) + cos(k2));
            d_y(ii,jj) = t*(sin(k1) + sin(k2));
            d_z(ii,jj) = 2*V - 4*lambda*( sin(k1) - sin(k2) - sin(k12));
            D = sqrt( d_x(ii,jj)^2 + d_y(ii,jj)^2 + d_z(ii,jj)^2 );
            %Hamiltonian + eigenvalue
            H = [d_z(ii,jj), d_x(ii,jj)-i*d_y(ii,jj); d_x(ii,jj)+i*d_y(ii,jj), -d_z(ii,jj)];
            [Psi,energy] = eig(H);
            U(:,:,ii,jj) = Psi; % 4-tensor
                Psi(1,:) = Psi(1,:)*exp(-i*0);
                Psi(2,:) = Psi(2,:)*exp(-i*(k1+k2)/2);
            E1(ii,jj) = energy(1,1);
            E2(ii,jj) = energy(2,2);
            if(k_y==0)
                E_kx(1,ii) = E1(ii,jj);
                E_kx(2,ii) = E2(ii,jj);
            end
            if(k_x==0)
                E_ky(1,jj) = E1(ii,jj);
                E_ky(2,jj) = E2(ii,jj);
            end
            %Normalized d-vector elements
            d_x(ii,jj) = d_x(ii,jj)/D;
            d_y(ii,jj) = d_y(ii,jj)/D;
            d_z(ii,jj) = d_z(ii,jj)/D;
        end
    end
%%
%Defining k-path (œcie¿ka k) & Energy along k-path
     grid = 200; %k_path grid
     parameters = [t, V, lambda]; %wektor z zadanymi parametrami
     k_path = K_gridGMKG(b1,b2,K,M,grid);
     Energy_k_path(k_path,parameters,a1,a2,grid);
%     Berry_k_path(k_path,parameters,grid,a1,a2);
        
     %Returns chern number for valence and conducting band
%     [C_valence,C_conducting] = Berry(U,kx,ky,parameters,a1,a2); 

%%
%Plotting D-vector
% figure(5);
% mesh(d_x,d_y,d_z);
% title(sprintf('D-vector \n using parameters: t =%1.1f, V = %1.2f and \\lambda = %1.2f ',t,V,lambda));
% xlabel(sprintf('d_x')); ylabel(sprintf('d_y')); zlabel(sprintf('d_z'));
% 
% %Plotting 3D plot E(kx,ky) and 2D plots (E(kx,0) & E(0,ky)
% figure(2); % 3D plot
% Colour = gradient(E1-E2);
% meshc(kx,ky,E1); hold on
% meshc(kx,ky,E2);
% xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
% xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% yticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
% yticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% xlabel(sprintf('k_x [\\pi/a]')); ylabel(sprintf('k_y [\\pi/a]')); zlabel('E [eV]');
% title(sprintf('Energy spectrum \n using parameters: t =%1.1f, V = %1.2f and \\lambda = %1.2f ',t,V,lambda))

% figure(3); %along x axis -> ky = 0
% plot(kx,E_kx(1,:),'k-', kx, E_kx(2,:),'k-');
% xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
% xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% xlabel(sprintf('k_x [\\pi/a]')); ylabel('E [eV]');
% title(sprintf('Energy spectrum along k_y = 0 \n using parameters: t =%1.1f, V = %1.2f and \\lambda = %1.2f ',t,V,lambda))
% 
% figure(4); % along y-axis  -> kx = 0
% plot(ky,E_ky(1,:),'k-', ky, E_ky(2,:),'k-');
% xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
% xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% xlabel(sprintf('k_y [\\pi/a]')); ylabel('E [eV]');
% title(sprintf('Energy spectrum along k_x = 0 \n using parameters: t =%1.1f, V = %1.2f and \\lambda = %1.2f ',t,V,lambda))
% 



