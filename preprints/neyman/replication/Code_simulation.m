%%% CODE FOR "A Neyman Orthogonalization Approach to the Incidental Parameter Problem in Likelihood Models", by Bonhomme, Jochmans, Weidner

%%% This code performs a Monte Carlo simulation
%%% It produces Tables 3 and 4 in the supplement

clear
clc
% Takes the team network from the empirical application
load('data_research_processed');

% Monte Carlo
tic
% Number of simulated samples
Nsim=300;

Results0=zeros(Nsim,4);
Results1=zeros(Nsim,4);
Results2=zeros(Nsim,4);
Results3=zeros(Nsim,4);
Results4=zeros(Nsim,4);
Results5=zeros(Nsim,4);
Results6=zeros(Nsim,4);

% parallel loop
parfor jsim=1:Nsim

    rng(10000+jsim)

    %%% DGP
    beta0=1;
    lambda0=log(beta0);
    gamma0=1;
    sigma01=1;
    sigma02=1;
    alpha0=exp(randn(size(A,2),1));
    logYmean=(sum(A,2)==2).*(lambda0+(gamma0^(-1))*log((1/2)*(A*(alpha0.^gamma0)))+sigma02*randn(size(A,1),1))...
        +(sum(A,2)==1).*(log(A*alpha0)+sigma01*randn(size(A,1),1));
    Ymean=exp(logYmean);

    y=Ymean;
    x=A;

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

    %%% Uncorrected estimator (q=0)

    delta_loop3=fminbnd(@(a) CES_lik_small(a,full(x2B),full(y2B),alpha_tilde), -500,500);


    Vect=log(y2B)-delta_loop3*log((1/2)*x2B*(alpha_tilde.^(1/delta_loop3)));

    lambda_loop3=mean(Vect);

    sigma2_loop3=mean((Vect-lambda_loop3).^2);

    par_loop3=[lambda_loop3 delta_loop3 log(sigma2_loop3)]';

    % variance in teams of size 1

    var_size1_loop3=mean((log(y1B)-log(x1B*alpha_tilde)).^2);

    par_hat=[1/par_loop3(2) par_loop3(1) exp(par_loop3(3)) var_size1_loop3];

    Results0(jsim,:)=par_hat;


    % prepare matrices - subnetworks of 2 workers and 3 teams of the form {A,B,A&B}

    Mat_alpha_prelim=zeros(0,2);
    Mat_meanoutput=zeros(0,3);

    for ii=1:size(x2B,1)

        zz=find(x2B(ii,:)==1);

        % preliminary estimates (in logs)

        Mat_alpha_prelim(ii,1)=log(alpha_tilde(zz(1)));

        Mat_alpha_prelim(ii,2)=log(alpha_tilde(zz(2)));

        % sole-authored article of worker 1

        vect1=y1B(x1B(:,zz(1))==1);
        Mat_meanoutput(ii,1)=log(vect1);

        % sole-authored article of worker 1

        vect2=y1B(x1B(:,zz(2))==1);
        Mat_meanoutput(ii,2)=log(vect2);

        % co-authored article of workers 1 and 2

        Mat_meanoutput(ii,3)=log(y2B(ii));

    end

    % minimization of projected score (q=1 to 6)

    warning('off')
    options=optimset('maxiter',5000,'Display','off');

    par_hat1=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,1),par_hat,options);

    par_hat2=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,2),par_hat1,options);

    par_hat3=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,3),par_hat2,options);

    par_hat4=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,4),par_hat3,options);

    par_hat5=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,5),par_hat4,options);

    par_hat6=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,6),par_hat5,options);


    % Store parameter estimates
    Results1(jsim,:)=par_hat1;
    Results2(jsim,:)=par_hat2;
    Results3(jsim,:)=par_hat3;
    Results4(jsim,:)=par_hat4;
    Results5(jsim,:)=par_hat5;
    Results6(jsim,:)=par_hat6;

end
toc

% Produce the numbers in Tables 3 and 4 in the supplement
disp([median(Results0(:,1)) mean(Results0(:,1)) quantile(Results0(:,1),0.025) quantile(Results0(:,1),0.975)...
    median(exp(Results0(:,2))) mean(exp(Results0(:,2))) quantile(exp(Results0(:,2)),0.025) quantile(exp(Results0(:,2)),0.975)...
    median(Results0(:,3)) mean(Results0(:,3)) quantile(Results0(:,3),0.025) quantile(Results0(:,3),0.975)...
    median(Results0(:,4)) mean(Results0(:,4)) quantile(Results0(:,4),0.025) quantile(Results0(:,4),0.975)])
disp([median(Results1(:,1)) mean(Results1(:,1)) quantile(Results1(:,1),0.025) quantile(Results1(:,1),0.975)...
    median(exp(Results1(:,2))) mean(exp(Results1(:,2))) quantile(exp(Results1(:,2)),0.025) quantile(exp(Results1(:,2)),0.975)...
    median(Results1(:,3)) mean(Results1(:,3)) quantile(Results1(:,3),0.025) quantile(Results1(:,3),0.975)...
    median(Results1(:,4)) mean(Results1(:,4)) quantile(Results1(:,4),0.025) quantile(Results1(:,4),0.975)])
disp([median(Results2(:,1)) mean(Results2(:,1)) quantile(Results2(:,1),0.025) quantile(Results2(:,1),0.975)...
    median(exp(Results2(:,2))) mean(exp(Results2(:,2))) quantile(exp(Results2(:,2)),0.025) quantile(exp(Results2(:,2)),0.975)...
    median(Results2(:,3)) mean(Results2(:,3)) quantile(Results2(:,3),0.025) quantile(Results2(:,3),0.975)...
    median(Results2(:,4)) mean(Results2(:,4)) quantile(Results2(:,4),0.025) quantile(Results2(:,4),0.975)])
disp([median(Results3(:,1)) mean(Results3(:,1)) quantile(Results3(:,1),0.025) quantile(Results3(:,1),0.975)...
    median(exp(Results3(:,2))) mean(exp(Results3(:,2))) quantile(exp(Results3(:,2)),0.025) quantile(exp(Results3(:,2)),0.975)...
    median(Results3(:,3)) mean(Results3(:,3)) quantile(Results3(:,3),0.025) quantile(Results3(:,3),0.975)...
    median(Results3(:,4)) mean(Results3(:,4)) quantile(Results3(:,4),0.025) quantile(Results3(:,4),0.975)])
disp([median(Results4(:,1)) mean(Results4(:,1)) quantile(Results4(:,1),0.025) quantile(Results4(:,1),0.975)...
    median(exp(Results4(:,2))) mean(exp(Results4(:,2))) quantile(exp(Results4(:,2)),0.025) quantile(exp(Results4(:,2)),0.975)...
    median(Results4(:,3)) mean(Results4(:,3)) quantile(Results4(:,3),0.025) quantile(Results4(:,3),0.975)...
    median(Results4(:,4)) mean(Results4(:,4)) quantile(Results4(:,4),0.025) quantile(Results4(:,4),0.975)])
disp([median(Results5(:,1)) mean(Results5(:,1)) quantile(Results5(:,1),0.025) quantile(Results5(:,1),0.975)...
    median(exp(Results5(:,2))) mean(exp(Results5(:,2))) quantile(exp(Results5(:,2)),0.025) quantile(exp(Results5(:,2)),0.975)...
    median(Results5(:,3)) mean(Results5(:,3)) quantile(Results5(:,3),0.025) quantile(Results5(:,3),0.975)...
    median(Results5(:,4)) mean(Results5(:,4)) quantile(Results5(:,4),0.025) quantile(Results5(:,4),0.975)])
disp([median(Results6(:,1)) mean(Results6(:,1)) quantile(Results6(:,1),0.025) quantile(Results6(:,1),0.975)...
    median(exp(Results6(:,2))) mean(exp(Results6(:,2))) quantile(exp(Results6(:,2)),0.025) quantile(exp(Results6(:,2)),0.975)...
    median(Results6(:,3)) mean(Results6(:,3)) quantile(Results6(:,3),0.025) quantile(Results6(:,3),0.975)...
    median(Results6(:,4)) mean(Results6(:,4)) quantile(Results6(:,4),0.025) quantile(Results6(:,4),0.975)])


% save results
save('Results_sim')