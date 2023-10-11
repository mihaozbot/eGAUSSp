function model = merge_rules(model, kappa_V_th)

%Unpack model
c = model.c;
mu = model.mu;
n = model.n;
S = model.S;
n_z = size(mu,1);

%%Cluster merging mechanism
merge = 1; %Merge until no more merging is required
while merge
    merge = 0;

    V = NaN(c,c);
    Sigma_ij = NaN(n_z,n_z,c,c);
    mu_ij = NaN(n_z,c,c);
    n_ij = NaN(c,c);

    %Compute volumes
    for i = 1:c

        %Compute self volumes
        V(i,i) = det(S(:,:,i)/n(i));

        for j = 1:c
            if i==j
                continue
            end


            %Combine i and j 
            n_ij(i,j) = n(i) + n(j);
            mu_ij(:,i,j) = (n(i)*mu(:,i) + n(j)*mu(:,j))/n_ij(i,j);
            ZiTZi = (n(i)-1)*(1/n(i))*S(:,:,i) + diag(mu(:,i))'*...
                (ones(n(i),n_z)')*ones(n(i),n_z)*diag(mu(:,i));
            ZjTZj = (n(j)-1)*(1/n(j))*S(:,:,j) + diag(mu(:,j))'*...
                (ones(n(j),n_z)')*ones(n(j),n_z)*diag(mu(:,j));
            Sigma_ij(:,:,i,j) = (1/(n_ij(i,j)-1))*(ZiTZi + ZjTZj ...
                -diag(mu_ij(:,i,j))'*(ones(n_ij(i,j),n_z)')*...
                ones(n_ij(i,j),n_z)*diag(mu_ij(:,i,j)));
        
            %Compute volume
            V(i,j) = det(Sigma_ij(:,:,i,j));

            if  V(i,j) < 0
                V(i,j) = NaN; %Numerical instability safeguard
            end

        end
    end

    %Compute relation
    kappa_V = Inf(c,c);
    for i = 1:c
        for j = 1:c
            if i==j
                continue;
            end
            kappa_V(i,j) = V(i,j)/(V(i,i)+ V(j,j));
        end
    end

    %Check condition
    cond = (kappa_V < kappa_V_th) ;
    if ~any(any(cond))
        continue;
    end

    %Compute indexes
    merge_check = nan(c,c);
    merge_check(cond) = kappa_V(cond);
    [merge_check,i_e_min] = min(merge_check);
    [~, j_merge] = min(merge_check);
    i_merge = i_e_min(j_merge);

    %Merge i and j
    n(i_merge) = n_ij(i_merge,j_merge);
    mu(:,i_merge) = mu_ij(:,i_merge,j_merge);
    S(:,:,i_merge) = Sigma_ij(:,:,i_merge,j_merge)*n_ij(i_merge,j_merge);

    %...and remove j
    n(j_merge) = [];
    mu(:,j_merge) = [];
    S(:,:,j_merge) = [];
    c = c - 1;
    merge = 1;

end

%Pach model
model.c = c;
model.mu = mu;
model.n = n;
model.S = S;

end

