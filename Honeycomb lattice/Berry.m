function [C_valence,C_conducting] = Berry(U,kx,ky,parameters,a1,a2)
%Calculating berry phase
t = parameters(1); V = parameters(2); lambda = parameters(3);
F_C = zeros(length(kx),length(ky));
F_V = zeros(length(kx),length(ky));
Chern1 = 0; Chern2 = 0; p=1;
KX = 4*pi/9*sqrt(3);
Kax = [KX,KX/2,-KX/2,-KX,-KX/2,KX/2];
Kay = [0,2*pi/3,2*pi/3,0,-2*pi/3,-2*pi/3];
for jj=1:length(ky)-1
    for ii=1:length(kx)-1
        x12 = dot(U(:,1,ii,jj),U(:,1,ii+1,jj));
        x23 = dot(U(:,1,ii+1,jj),U(:,1,ii+1,jj+1));
        x34 = dot(U(:,1,ii+1,jj+1),U(:,1,ii,jj+1));
        x41 = dot(U(:,1,ii,jj+1),U(:,1,ii,jj));
        berry_C = angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
        
        x12 = dot(U(:,2,ii,jj),U(:,2,ii+1,jj));
        x23 = dot(U(:,2,ii+1,jj),U(:,2,ii+1,jj+1));
        x34 = dot(U(:,2,ii+1,jj+1),U(:,2,ii,jj+1));
        x41 = dot(U(:,2,ii,jj+1),U(:,2,ii,jj));
        berry_V = angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
        
        d2k = ( (kx(ii)-kx(ii+1))^2 + (ky(jj)-ky(jj))^2 )^(1/2)*( (kx(ii+1)-kx(ii))^2 + (ky(jj+1)-ky(jj+1))^2 )^(1/2);
        
        F_C(ii,jj) = berry_C/d2k; 
        F_V(ii,jj) = berry_V/d2k;
 
        %Calculating the 1st Chern number for both bands
        if((abs(kx(ii))<4*sqrt(3)*pi/9)&&(abs(ky(jj))<2*pi/3))
            if((ky(jj)<(sqrt(3)*kx(ii)+4/3*pi))&&(ky(jj)>(sqrt(3)*kx(ii)-4/3*pi)))
              if((ky(jj)<(-sqrt(3)*kx(ii)+4/3*pi))&&(ky(jj)>(-sqrt(3)*kx(ii)-4/3*pi)))
%                  if(inpolygon(kx,ky,Kax,Kay)==true)
                    Chern1=Chern1+berry_C;
                    Chern2=Chern2+berry_V;
%                  end   
              end
            end    
        end
        
    end
end
%%
%Berry phase&curvature plot
C_valence = Chern1/2/pi
C_conducting = Chern2/2/pi


figure(8);
meshc(kx,ky,F_C);   
title(sprintf('Berry curvature for conducting band \n using parameters: V =%1.2f, t = %1.1f and \\lambda = %1.2f ',V,t,lambda));
xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
yticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
yticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
xlabel(sprintf('k_x [\\pi/a]')); ylabel(sprintf('k_y [\\pi/a]'));

figure(10);
meshc(kx,ky,F_V);
title(sprintf('Berry curvature for valence band \n using parameters: V =%1.2f, t = %1.1f and \\lambda = %1.2f ',V,t,lambda));
xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
yticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
yticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
xlabel(sprintf('k_x [\\pi/a]')); ylabel(sprintf('k_y [\\pi/a]'));

end

