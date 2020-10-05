function Berry_kagome(U,kx,ky,parameters,a1,a2)
%Calculating berry phase
t1 = parameters(1); L1 = parameters(2); 
t2 = parameters(3); L2 = parameters(4);
    KX = 4*pi/9*sqrt(3);
    Kax = [KX,KX/2,-KX/2,-KX,-KX/2,KX/2,KX]/norm(a1);
    Kay = [0,2*pi/3,2*pi/3,0,-2*pi/3,-2*pi/3,0]/norm(a2);
    F_1 = zeros(length(kx),length(ky));
    F_2 = zeros(length(kx),length(ky));
    F_3 = zeros(length(kx),length(ky));
    Chern1 = 0; Chern2 = 0; Chern3 = 0;
    for jj=1:length(ky)-1
        for ii=1:length(kx)-1
                    x =    dot(U(:,1,ii,jj),U(:,1,ii+1,jj))/abs(dot(U(:,1,ii,jj),U(:,1,ii+1,jj))) ;
                    x = x *dot(U(:,1,ii+1,jj),U(:,1,ii+1,jj+1))/abs(dot(U(:,1,ii+1,jj),U(:,1,ii+1,jj+1))) ;
                    x = x *dot(U(:,1,ii+1,jj+1),U(:,1,ii,jj+1))/abs(dot(U(:,1,ii+1,jj+1),U(:,1,ii,jj+1)));
                    x = x *dot(U(:,1,ii,jj+1),U(:,1,ii,jj))/abs(dot(U(:,1,ii,jj+1),U(:,1,ii,jj)));
                    berry_1 = angle(x);

                    y = dot(U(:,2,ii,jj),U(:,2,ii+1,jj))/abs(dot(U(:,2,ii,jj),U(:,2,ii+1,jj)));
                    y = y*dot(U(:,2,ii+1,jj),U(:,2,ii+1,jj+1))/abs(dot(U(:,2,ii+1,jj),U(:,2,ii+1,jj+1)));
                    y = y*dot(U(:,2,ii+1,jj+1),U(:,2,ii,jj+1))/abs(dot(U(:,2,ii+1,jj+1),U(:,2,ii,jj+1)));
                    y = y*dot(U(:,2,ii,jj+1),U(:,2,ii,jj))/abs(dot(U(:,2,ii,jj+1),U(:,2,ii,jj)));
                    berry_2 = angle(y);

                    y = dot(U(:,3,ii,jj),U(:,3,ii+1,jj))/abs(dot(U(:,3,ii,jj),U(:,3,ii+1,jj)));
                    y = y*dot(U(:,3,ii+1,jj),U(:,3,ii+1,jj+1))/abs(dot(U(:,3,ii+1,jj),U(:,3,ii+1,jj+1)));
                    y = y*dot(U(:,3,ii+1,jj+1),U(:,3,ii,jj+1))/abs(dot(U(:,3,ii+1,jj+1),U(:,3,ii,jj+1)));
                    y = y*dot(U(:,3,ii,jj+1),U(:,3,ii,jj))/abs(dot(U(:,3,ii,jj+1),U(:,3,ii,jj)));
                    berry_3 = angle(y);

                    d2k = ( (kx(ii)-kx(ii+1))^2 + (ky(jj)-ky(jj))^2 )^(1/2)*( (kx(ii+1)-kx(ii))^2 + (ky(jj+1)-ky(jj+1))^2 )^(1/2);

                    F_1(ii,jj) = berry_1/d2k; 
                    F_2(ii,jj) = berry_2/d2k; 
                    F_3(ii,jj) = berry_3/d2k;

                    %Calculating the 1st Chern number for both bands
                    if((abs(kx(ii))<4*pi/3&&(abs(ky(jj))<2*pi/3*sqrt(3))))
                        if((ky(jj)<(-sqrt(3)*kx(ii)+4*pi/sqrt(3)))&&(ky(jj)>(sqrt(3)*kx(ii)-4*pi)))
                          if((ky(jj)<(sqrt(3)*kx(ii)+4*pi/sqrt(3)))&&(ky(jj)>(-sqrt(3)*kx(ii)-4*pi/sqrt(3))))
                            %if(inpolygon(kx(ii),ky(jj),Kax,Kay)==true)
                              Chern1 = Chern1 + berry_1/2/pi;
                              Chern2 = Chern2 + berry_2/2/pi;
                              Chern3 = Chern3 + berry_3/2/pi;
                            %end
                          end
                        end    
                    end

        end
    end

    %%
    C1 = Chern1
    C2 = Chern2
    C3 = Chern3
%Berry phase&curvature plot
% figure(8);
% meshc(kx,ky,F_1); hold on
% plot(Kax,Kay,'k-','Linewidth',2);hold off
% title(sprintf('Berry curvature for 1st band \n using parameters: t1 =%1.1f, t2 = %1.2f, \\lambda1 = %1.2f and \\lambda2 = %1.2f ',t1,t2,L1,L2));
% xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
% xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% yticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
% yticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% xlabel(sprintf('k_x [\\pi/a]')); ylabel(sprintf('k_y [\\pi/a]'));
% 
% figure(10);
% meshc(kx,ky,F_2);hold on
% plot(Kax,Kay,'k-','Linewidth',2);hold off
% title(sprintf('Berry curvature for 2nd band \n using parameters: t1 =%1.1f, t2 = %1.2f, \\lambda1 = %1.2f and \\lambda2 = %1.2f ',t1,t2,L1,L2));
% xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
% xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% yticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
% yticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% xlabel(sprintf('k_x [\\pi/a]')); ylabel(sprintf('k_y [\\pi/a]'));
% 
% figure(17);
% meshc(kx,ky,F_3);hold on
% plot(Kax,Kay,'k-','Linewidth',2);hold off
% title(sprintf('Berry curvature for 3rd band \n using parameters: t1 =%1.1f, t2 = %1.2f, \\lambda1 = %1.2f and \\lambda2 = %1.2f ',t1,t2,L1,L2));
% xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
% xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% yticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
% yticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
% xlabel(sprintf('k_x [\\pi/a]')); ylabel(sprintf('k_y [\\pi/a]'));

end

