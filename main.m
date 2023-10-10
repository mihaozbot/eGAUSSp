%Miha Ožbot 2022

clc; clear all; close all;

times_sigma = 2; %Number of standard deviations that are plotted
N_plot = 1; %Number of samples in cluster before it is plotted
online_display = false; %Enable/Disable 1/0 display of clusters online
load_data = false; %Data used in the paper
dispaly_top_three = false; % Display just the top three clusters


%Measurements
if load_data

    N = 400;
    load('Z2D.mat');
    z = Z(:,1:N);
    
else

    z_1(1,:) = 1.*randn(1,0);
    z_1(2,:) = z_1(1,:) + 0.8.*randn(1,0);
    z_2(1,:) = 1.*randn(1,100);
    z_2(2,:)  = -z_2(1,:) + (8 + 0.85.*randn(1,100));
    z_3(1,:) = -1.5 + 0.5.*randn(1,150);
    z_3(2,:)  = -z_3(1,:) + (4 + 0.5.*randn(1,150));
    z = [z_1,z_2,z_3];
    N = size(z,2);
    z = z(:,randi([1,N],1,N));

end

%Initialization
m = size(z,1);
eGAUSSp.c = 1;
eGAUSSp.mu = z(:,1);
eGAUSSp.n(1) = 1;
eGAUSSp.S(:,:,1) = zeros(m,m,1);

%Parameters
par.kappa_join = 1.5;
N_r = 4;
d_max = zeros(m,1);
for i = 1:m
    d_max(i) = (max(z(i,:))-min(z(i,:)))/(2*N_r);
end
par.Gamma_max = exp(-3^2); exp(-min(d_max)^2);
par.N_max = N_r;
par.S_0 = 1e-2;

for k = 2:N

    eGAUSSp = evolve_egaussp(z(:,k),eGAUSSp,par);

    %Online clustering dispaly
    if online_display
        display_clusters
    end

end

if dispaly_top_three
    eGAUSSp.c = 3;
    [~,idx] = maxk(eGAUSSp.n,eGAUSSp.c);
    eGAUSSp.n = eGAUSSp.n(idx);
    eGAUSSp.mu = eGAUSSp.mu(:,idx);
    eGAUSSp.S = eGAUSSp.S(:,:,idx);
end

save('eGAUSSp')

%Final display
display_clusters
