function eGAUSSp = evolve_egaussp(z,eGAUSSp,par)
    
    %Unpack model
    c = eGAUSSp.c;
    n = eGAUSSp.n;
    mu = eGAUSSp.mu;
    S = eGAUSSp.S;
    m = size(mu,1);
    
    %Unpack method constants
    N_max = par.N_max;
    Gamma_max = par.Gamma_max;
    kappa_join= par.kappa_join;
    S_0 = par.S_0;

    %Compute distance
    Gamma = zeros(c,1);
    d2 = zeros(c,1);

    for i = 1:c
        if n < N_max
             d2(i) = ((z- mu(:,i))'*(z - mu(:,i))); %Euclidian distance
        else
            d2(i) = ((z- mu(:,i))'*pinv(S(:,:,i)/n(i))*(z - mu(:,i))); %Mahalanobis distance
        end
        Gamma(i) = exp(-d2(i));

    end

    %Computer maximum activation
    [~,j] = max(Gamma);

    if Gamma(j) > Gamma_max

        %Increment cluster
        e = z - mu(:,j); %Distance of new data from center
        mu(:,j) = mu(:,j)+ 1/(1 + n(j))*e; %Center update
        S(:,:,j)  = S(:,:,j) + e*(z - mu(:,j))'; %un-normalized covariance matrix
        n(j) = n(j) + 1; %Increase number of samples in cluster

    else

        %Add and initialize new cluster
        c = c + 1;
        n(c) = 1;
        mu(:,c) = z;
        S(:,:,c) = S_0*eye(m);
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

    %Pack the eGAUSS+ model
    eGAUSSp.c = c;
    eGAUSSp.n = n;
    eGAUSSp.mu = mu;
    eGAUSSp.S = S;
