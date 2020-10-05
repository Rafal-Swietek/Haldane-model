function Energy_k_path_kagome(k_path,parameters,a1,a2)
%Function solving the eigenvalue problem along the Gamma-K-M-Gamma path and
%   returning a matrix of eigenvecotrs along this path
%U_k_path = zeros(2,2,length(k_path)); %matrix for eigenvectors (array of matrices of eigenvectors for different k)
%Defining k-path (œcie¿ka k) & Energy along k-path
    t1 = parameters(1);
    L1 = parameters(2);
    t2 = parameters(3);
    L2 = parameters(4);
    E1 = zeros(1,length(k_path));
    E2 = zeros(1,length(k_path)); % last argument defines the discretization of the track Gamma-K-M-Gamma
    E3 = zeros(1,length(k_path));
    for ii=1:length(k_path)
            k_x = k_path(1,ii); k_y = k_path(2,ii);
            k1 = k_x*a1(1)+k_y*a1(2);    k2 = k_x*a2(1)+k_y*a2(2);
            u = exp(i*k1);      v = exp(i*k2); %conj(u) = 1/u, because its exponential
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
            energy = eig(H);
            %U(:,:,ii,jj) = Psi; % 4-tensor
            E1(ii) = energy(1);
            E2(ii) = energy(2); 
            E3(ii) = energy(3);
    end
    
% %   Analitycal eigenvectors along k-path
%     fi =zeros(1,length(k_path)); theta = zeros(1,length(k_path)); 
%     [fi,theta,d] = cart2sph(d_x,d_y,d_z);
%     for ii=1:length(k_path) %analitycal eigenvectors
%           U_k_path(1,1,ii) = -exp(-i*fi(ii))*cos(theta(ii)/2);
%           U_k_path(2,1,ii) = -sin(theta(ii)/2);
%           U_k_path(1,2,ii) = -exp(-i*fi(ii))*sin(theta(ii)/2);
%           U_k_path(2,2,ii) = cos(theta(ii)/2);
%     end
    figure(67);
    %Plotting energy solution
    y = -15:1:15; x1 = length(k_path)/2.5*ones(1,length(y))+1; x2 = length(k_path)/2.5*1.5*ones(1,length(y))+1; %dashed vertical line
    %y = min(E1):1:max(E2); x1 = length(k_path)/3*ones(1,length(y)); x2 = length(k_path)/3*2*ones(1,length(y)); %dashed vertical line
    plot(x1,y,'k--',x2,y,'k--',1:length(k_path),E1,'k-', 1:length(k_path),E2,'k-',1:length(k_path),E3,'k-'); drawnow
    xticks([0 length(k_path)/2.5+1 1.5*length(k_path)/2.5+1 length(k_path)]);
    %xticks([0 length(k_path)/3 2*length(k_path)/3 length(k_path)]);
    xticklabels({'\Gamma','K','M','\Gamma'});
    axis([0 length(k_path) min(E1) max(E3)]);
    ylabel('E [eV]');
    title(sprintf('Energy spectrum for given k-path:\\Gamma-K-M-\\Gamma \n using parameters: t1 =%1.1f, t2 = %1.2f, \\lambda1 = %1.2f and \\lambda2 = %1.2f ',t1,t2,L1,L2));
%     
%    % Creating .gif file for animation
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind, cm] = rgb2ind(im,256);
%     if L == -L_0;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
end

