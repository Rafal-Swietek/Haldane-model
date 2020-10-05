function Energy_k_path(k_path,parameters,a1,a2,grid)
%Function solving the eigenvalue problem along the Gamma-K-M-Gamma path and
%   returning a matrix of eigenvecotrs along this path
%U_k_path = zeros(2,2,length(k_path)); %matrix for eigenvectors (array of matrices of eigenvectors for different k)
%Defining k-path (œcie¿ka k) & Energy along k-path
%L_0 = 1;
%figure(1);
    t = parameters(1);
    V = parameters(2);
    lambda = parameters(3);
    E1 = zeros(1,length(k_path));
    E2 = zeros(1,length(k_path)); % last argument defines the discretization of the track Gamma-K-M-Gamma
    for ii=1:length(k_path)
       kx = k_path(1,ii);
       ky = k_path(2,ii);
       d_x = t*(1 + cos(kx*a1(1)+ky*a1(2)) + cos(kx*a2(1)+ky*a2(2)));
       d_y = t*( sin(kx*a1(1)+ky*a1(2)) + sin(kx*a2(1)+ky*a2(2)) );
       d_z = 2*V - 4*lambda*( sin(kx*a1(1)+ky*a1(2)) - sin(kx*a2(1)+ky*a2(2)) - sin(kx*(a1(1)-a2(2))+ky*(a1(2)-a2(2))) );
       H = [d_z, d_x-i*d_y; d_x+i*d_y, -d_z];
       energy = eig(H);
%        [Psi,energy] = eig(H);
%        U_k_path(:,:,ii) = Psi;
       E1(ii) = energy(1);
       E2(ii) = energy(2);
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
    
    %Plotting energy solution
    y = -15:1:15; x1 = length(k_path)/2.5*ones(1,length(y))+1; x2 = length(k_path)/2.5*1.5*ones(1,length(y))+1; %dashed vertical line
    %y = min(E1):1:max(E2); x1 = length(k_path)/3*ones(1,length(y)); x2 = length(k_path)/3*2*ones(1,length(y)); %dashed vertical line
    plot(x1,y,'k--',x2,y,'k--',1:length(k_path),E1,'k-', 1:length(k_path),E2,'k-');
    xticks([0 length(k_path)/2.5+1 1.5*length(k_path)/2.5+1 length(k_path)]);
    %xticks([0 length(k_path)/3 2*length(k_path)/3 length(k_path)]);
    xticklabels({'\Gamma','K','M','\Gamma'});
    axis([0 length(k_path) min(E1) max(E2)]); hold on
    ylabel('E [eV]');
    title(sprintf('Energy spectrum for given k-path:\\Gamma-K-M-\\Gamma \n using parameters: t =%1.1f , V=%1.2f and \\lambda = %1.2f ',t,V,lambda));
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

