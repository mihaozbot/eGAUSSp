c = eGAUSSp.c;
n = eGAUSSp.n;
mu = eGAUSSp.mu;
S = eGAUSSp.S;

angle = 0:pi/100:2*pi; %Angles around a circle
xy = zeros(length(angle),2,c); %Points on the ellipse
Sigma = zeros(2,2,c); %Normalized covariance matrix

for i = 1:1:c
    Sigma(:,:,i) = S(:,:,i)/n(i);
    [eigvec, eigval] = eig(Sigma(:,:,i));
    xy(:,:,i) = [cos(angle'),sin(angle')] *times_sigma*sqrt(eigval) * eigvec';
end

color = lines(c);

figure(1); hold off;
plot(z(1,1:k),z(2,1:k),'r.','markersize',10); hold on;
xlabel('z_1')
ylabel('z_2')
title('eGAUSS+')

for i = 1:1:c
    if (n(i) >= N_plot)
        plot(mu(1,i),mu(2,i),'o','Color','k','markersize',2,'linewidth',3)
        if (det(Sigma(:,:,i)) > 1e-15)
            plot(xy(:,1,i) + mu(1,i),xy(:,2,i)+ mu(2,i),'Color','b','linewidth',1.5)
        end
    end
end

