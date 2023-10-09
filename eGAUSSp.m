%Miha Ožbot 2022

clc; clear all; close all;

times_sigma = 2; %Number of standard deviations that are plotted
N_plot = 1; %Number of samples in cluster before it is plotted
online_display = false; %Enable/Disable 1/0 display of clusters online
load_data = false; %Data used in the paper

%Measurements
if load_data

    N = 400;
    load('Z2D.mat');
    z = Z(:,1:N);
    
else

    z_1(1,:) = 1.*randn(1,250);
    z_1(2,:) = z_1(1,:) + 0.8.*randn(1,250);
    z_2(1,:) = 1.*randn(1,100);
    z_2(2,:)  = -z_2(1,:) + (8 + 0.85.*randn(1,100));
    z_3(1,:) = -1.5 + 0.5.*randn(1,50);
    z_3(2,:)  = -z_3(1,:) + (4 + 0.5.*randn(1,50));
    z = [z_1,z_2,z_3];
    N = size(z,2);
    z = z(:,randi([1,N],1,N));

end

%Initialization
c = 1;
m = size(z,1);
mu = z(:,1);
n(1) = 1;
S(:,:,1) = zeros(m,m,1);

%Parameters
kappa_join = 1.4;
N_r = 4;
d_max = zeros(m,1);
for i = 1:m
    d_max(i) = (max(z(i,:))-min(z(i,:)))/(2*N_r);
end
Gamma_max = exp(-min(d_max)^2);
N_max = 4;

for k = 2:N

    %Compute distance
    Gamma = zeros(c,1);
    d2 = zeros(c,1);

    for i = 1:c

        if n(i) < N_max
            d2(i) = ((z(:,k) - mu(:,i))'*(z(:,k) - mu(:,i))); %Euclidean distance
        else
            d2(i) = ((z(:,k)- mu(:,i))'*pinv(S(:,:,i)/n(i))*(z(:,k) - mu(:,i))); %Mahalanobis distance
        end

        Gamma(i) = exp(-d2(i));

    end

    %Computer maximum activation
    [~,j] = max(Gamma);

    if Gamma(j) > Gamma_max

        %Increment cluster
        e = z(:,k) - mu(:,j); %Distance of new data from center
        mu(:,j) = mu(:,j)+ 1/(1 + n(j))*e; %Center update
        S(:,:,j)  = S(:,:,j) + e*(z(:,k) - mu(:,j))'; %un-normalized covariance matrix
        n(j) = n(j) + 1; %Increase number of samples in cluster

    else

        %Add and initialize new cluster
        c = c + 1;
        n(c) = 1;
        mu(:,c) = z(:,k);
        S(:,:,c) = zeros(m,m,1);
        Gamma(c) = 1;

    end

    %Cluster merging mechanism
    merge = 1; %Merge until no more merging is required
    while merge

        V = NaN(c,c);
        Sigma_ij = NaN(m,m,c,c);
        mu_ij = NaN(m,c,c);
        n_ij = NaN(c,c);

        for i = 1:c

            if Gamma(i) > Gamma_max/4

                V(i,i) = (2*pi^(m/2)/(m*gamma(m/2)))*prod(eig(S(:,:,i)/n(i)));

                for j = i+1:c

                    if Gamma(j) > Gamma_max/4

                        n_ij(i,j) = n(i) + n(j);

                        mu_ij(:,i,j) = (n(i)*mu(:,i) + n(j)*mu(:,j))/n_ij(i,j);
                        ZiTZi = (n(i)-1)*(1/n(i))*S(:,:,i) + diag(mu(:,i))'*...
                            (ones(n(i),m)')*ones(n(i),m)*diag(mu(:,i));

                        ZjTZj = (n(j)-1)*(1/n(j))*S(:,:,j) + diag(mu(:,j))'*...
                            (ones(n(j),m)')*ones(n(j),m)*diag(mu(:,j));

                        Sigma_ij(:,:,i,j) = (1/(n_ij(i,j)-1))*(ZiTZi + ZjTZj ...
                            -diag(mu_ij(:,i,j))'*(ones(n_ij(i,j),m)')*...
                            ones(n_ij(i,j),m)*diag(mu_ij(:,i,j)));

                        V(i,j) = (2*pi^(m/2)/m*Gamma(m/2))*prod(eig(Sigma_ij(:,:,i,j)));

                        if  V(i,j) < 0
                            V(i,j) = NaN; %Numerical instability safeguard
                        end

                    end
                end
            end
        end

        kappa = NaN(c,c);
        for i = 1:c
            for j = i+1:c
                kappa(i,j) = V(i,j)/(V(i,i)+ V(j,j));
            end
        end

        [kappa_min,i_kappa_min] = min(kappa,[],'all');
        [i,j] = ind2sub(size(kappa),i_kappa_min);

        if kappa_min < kappa_join

            n(i) = n_ij(i,j);
            mu(:,i) = mu_ij(:,i,j);
            S(:,:,i) = Sigma_ij(:,:,i,j)*n_ij(i,j);

            n(j) = [];
            mu(:,j) = [];
            S(:,:,j) = [];
            Gamma(j) = [];

            c = c - 1;
            merge = 1;

        else
            merge = 0;
        end

    end

    %Online clustering dispaly
    if online_display
        display_clusters
    end

end


c = 3;
[~,idx] = maxk(n,c);
n = n(idx);
mu = mu(:,idx);
S = S(:,:,idx);

save('eGAUSSp','n','mu','S')

%Final display
display_clusters



