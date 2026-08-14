%%% CODE FOR "A Neyman Orthogonalization Approach to the Incidental Parameter Problem in Likelihood Models", by Bonhomme, Jochmans, Weidner

%%% This code produces the GMM results in the supplement

clear
clc

load('data_research_processed');
Ymean_local=Ymean;

% If MC=1 this is a Monte Carlo
% If MC=0 this is estimation on the data

MC=1;

% Loop for cross-fitting (or simulations if MC=1)

if MC==1
    Ncf=300;
else
    Ncf=100;
end

Res=zeros(Ncf,1);

tic

parfevalOnAll(@() warning('off','MATLAB:nearlySingularMatrix'), 0);
parfor jcf=1:Ncf

    rng(300+jcf)

    % simulation design (if MC=1)

    if MC==1

        beta0=1;
        lambda0=log(beta0);
        gamma0=1;
        sigma01=1/5;
        sigma02=1/5;
        alpha0=exp(randn(size(A,2),1));
        logYmean=(sum(A,2)==2).*(lambda0+(gamma0^(-1))*log((1/2)*(A*(alpha0.^gamma0)))+sigma02*randn(size(A,1),1))...
            +(sum(A,2)==1).*(log(A*alpha0)+sigma01*randn(size(A,1),1));
        Ymean=exp(logYmean);

    end

    if MC==1
        y=Ymean;
        x=A;
    else
        y=Ymean_local;
        x=A;
    end


    % Select all publications in 1-teams but 1 per author (at random)

    x1=x(sum(x,2)==1,:);

    vect_sel1=[];
    vect_sel2=[];


    for ii=1:size(x1,2)
        zz=find(x1(:,ii));
        nzz=size(zz,1);
        uu=randperm(nzz);
        vect_sel1=[vect_sel1;zz(uu(1:nzz-1))];
        vect_sel2=[vect_sel2;zz(uu(nzz))];
    end

    x1_holdout=x1(vect_sel1,:);

    y1=y(sum(x,2)==1);

    y1_holdout=y1(vect_sel1);

    % preliminary estimates

    alpha_tilde=exp(sum(log(y1_holdout).*x1_holdout)'./sum(x1_holdout)');

    % Select out one observation in each 1-team

    x1_keep=x1(vect_sel2,:);
    x_keep=[x1_keep;x(sum(x,2)==2,:)];

    y1_keep=y1(vect_sel2);
    y_keep=[y1_keep;y(sum(x,2)==2)];

    x=x_keep;
    y=y_keep;

    % Select the 2-worker teams

    x2=x(sum(x,2)==2,:);

    y2=y(sum(x,2)==2);

    % Select the 1-worker teams (for later)

    x1=x(sum(x,2)==1,:);

    y1=y(sum(x,2)==1);

    % keep only sole-authored papers of workers who co-author with others
    % we will select subnetworks of 2 workers and 3 teams of the form {A,B,A&B}

    vec_sub1=[];


    for ii=1:size(x2,1)

        zz=find(x2(ii,:)==1);

        ind1=find(x1(:,zz(1))==1);

        ind2=find(x1(:,zz(2))==1);

        vec_sub1=[vec_sub1;ind1;ind2];

    end

    vec_sub1=unique(vec_sub1);

    x2B=x2;

    y2B=y2;

    x1B=x1(vec_sub1,:);

    y1B=y1(vec_sub1,1);


    % construction of instruments for GMM
    KK=1;
    inst=zeros(size(x2B,1),KK);

    for ii=1:size(x2B,1)
        zz=find(x2B(ii,:)==1);
        for kk=1:KK
            inst(ii,kk)=(alpha_tilde(zz(1)).*alpha_tilde(zz(2))).^(kk/KK);
        end
    end

    Minst = inst - x2*(x2\inst);

    grid_gamma=(-2:.001:3)';
    res_obj=zeros(size(grid_gamma,1),1);
    for j_gamma=1:size(grid_gamma,1)
        res_obj(j_gamma)=sum((Minst'*(y2.^grid_gamma(j_gamma))).^2);
    end

    % minimization
    TF = islocalmin(res_obj);

    if sum(TF)~=2
        % case without two local minima
        % there is always a local minimum at zero
        Res(jcf)=-1000;
    else
        mm=find(TF==1);
        if mm(1)==2001
            % this case corresponds to the local arg minimum = 0
            Res(jcf)=grid_gamma(mm(2));
        else
            Res(jcf)=grid_gamma(mm(1));
        end
    end

end

toc

% mean
disp(mean(Res(Res~=-1000)))

% median
disp(median(Res(Res~=-1000)))

% standard deviation
disp(std(Res(Res~=-1000)))

% occurances with ony one local minimum (at 0)
disp(sum(Res==-1000))