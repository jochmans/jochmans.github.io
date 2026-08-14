%%% CODE FOR "A Neyman Orthogonalization Approach to the Incidental Parameter Problem in Likelihood Models", by Bonhomme, Jochmans, Weidner

%%% This code produces the diagnostic for q in Table 1 in the paper

clear
clc

tic

% Test statistic in the original sample
load('Results_est');
Allpar=zeros(6,4);
Allpar(1,:)=mean([Results_par1(:,1) exp(Results_par1(:,2)) Results_par1(:,3) Results_par1(:,4)]);
Allpar(2,:)=mean([Results_par2(:,1) exp(Results_par2(:,2)) Results_par2(:,3) Results_par2(:,4)]);
Allpar(3,:)=mean([Results_par3(:,1) exp(Results_par3(:,2)) Results_par3(:,3) Results_par3(:,4)]);
Allpar(4,:)=mean([Results_par4(:,1) exp(Results_par4(:,2)) Results_par4(:,3) Results_par4(:,4)]);
Allpar(5,:)=mean([Results_par5(:,1) exp(Results_par5(:,2)) Results_par5(:,3) Results_par5(:,4)]);
Allpar(6,:)=mean([Results_par6(:,1) exp(Results_par6(:,2)) Results_par6(:,3) Results_par6(:,4)]);

load('data_research_processed');

% Cross fitting
Ncf=10;

Res_Score_q=zeros(Ncf,4,5);
Res_Score_q1=zeros(Ncf,4,5);


parfevalOnAll(@() warning('off','MATLAB:nearlySingularMatrix'), 0);
parfor jcf=1:Ncf

    rng(1000+jcf)

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

        % sole-authored article of worker 2

        vect2=y1B(x1B(:,zz(2))==1);
        Mat_meanoutput(ii,2)=log(vect2);

        % co-authored article of workers 1 and 2

        Mat_meanoutput(ii,3)=log(y2B(ii));

    end

    % loop on orthogonalization order
    for qq=1:5
        % q+1-orthogonal estimate
        par0=Allpar(qq+1,:);
        par0(2)=log(par0(2));
        % score u_q
        Score_q=compute_score(par0,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,qq)/sqrt(size(Mat_meanoutput,1));
        % score u_q+1
        Score_q1=compute_score(par0,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,qq+1)/sqrt(size(Mat_meanoutput,1));
        Res_Score_q(jcf,:,qq)=Score_q;
        Res_Score_q1(jcf,:,qq)=Score_q1;
    end
end

% record score difference
Test_stat = sum(Res_Score_q-Res_Score_q1);
Test_stat_mat = squeeze(Test_stat);
Test_stat_mat = Test_stat_mat.';








load('data_research_processed');

Nboot=200;

Results_boot=zeros(Nboot,15);
Results_0_boot=zeros(Nboot,2);
Results_1_boot=zeros(Nboot,2);
Results_2_boot=zeros(Nboot,2);
Results_3_boot=zeros(Nboot,2);
Results_4_boot=zeros(Nboot,2);
Results_5_boot=zeros(Nboot,2);
Results_6_boot=zeros(Nboot,2);
Results_par0_boot=zeros(Nboot,4);
Results_par1_boot=zeros(Nboot,4);
Results_par2_boot=zeros(Nboot,4);
Results_par3_boot=zeros(Nboot,4);
Results_par4_boot=zeros(Nboot,4);
Results_par5_boot=zeros(Nboot,4);
Results_par6_boot=zeros(Nboot,4);

% Loop for parametric bootstrap
% We simulate DGPs under 6-order orthogonalized estimates, with the
% eta_hat's estimated based on sole-authored publications

Results_Test_stat_mat_boot_par1=zeros(Nboot,5);
Results_Test_stat_mat_boot_par2=zeros(Nboot,5);
Results_Test_stat_mat_boot_par3=zeros(Nboot,5);
Results_Test_stat_mat_boot_par4=zeros(Nboot,5);


parfor jboot=1:Nboot

    rng(200+jboot)

    % parametric simulation
    gamma0=Allpar(6,1);
    beta0=Allpar(6,2);
    lambda0=log(beta0);
    sigma02=sqrt(Allpar(6,3));
    sigma01=sqrt(Allpar(6,4));

    y1=Ymean(sum(A,2)==1);
    x1=A(sum(A,2)==1,:);

    alpha0=exp(sum(log(y1).*x1)'./sum(x1)');


    logYmean=(sum(A,2)==2).*(lambda0+(gamma0^(-1))*log((1/2)*(A*(alpha0.^gamma0)))+sigma02*randn(size(A,1),1))...
        +(sum(A,2)==1).*(log(A*alpha0)+sigma01*randn(size(A,1),1));
    y0=exp(logYmean);
    x0=A;



    % Loop for cross-fitting

    Ncf=10;

    Res_Score_q_boot=zeros(Ncf,4,5);
    Res_Score_q1_boot=zeros(Ncf,4,5);



    for jcf=1:Ncf



        y=y0;
        x=x0;


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

        for qq=1:5
            % q+1-orthogonal estimates
            par0=Allpar(qq+1,:);
            par0(2)=log(par0(2));
            % score u_q
            Score_q=compute_score(par0,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,qq)/sqrt(size(Mat_meanoutput,1));
            % score u_q+1
            Score_q1=compute_score(par0,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,qq+1)/sqrt(size(Mat_meanoutput,1));
            Res_Score_q_boot(jcf,:,qq)=Score_q;
            Res_Score_q1_boot(jcf,:,qq)=Score_q1;
        end

    end

    % Store cross-fitted estimates

    % Record bootstrapped score difference
    Test_stat_boot = sum(Res_Score_q_boot-Res_Score_q1_boot);
    Test_stat_mat_boot = squeeze(Test_stat_boot);
    Test_stat_mat_boot = Test_stat_mat_boot.';

    Results_Test_stat_mat_boot_par1(jboot,:)=Test_stat_mat_boot(:,1)';
    Results_Test_stat_mat_boot_par2(jboot,:)=Test_stat_mat_boot(:,2)';
    Results_Test_stat_mat_boot_par3(jboot,:)=Test_stat_mat_boot(:,3)';
    Results_Test_stat_mat_boot_par4(jboot,:)=Test_stat_mat_boot(:,4)';


end

toc


%%% LAST COLUMN OF TABLE 1
% Test statistic based on inverse variance matrix

% loop on q

for qq=1:5
    mat=[Results_Test_stat_mat_boot_par1(:,qq) Results_Test_stat_mat_boot_par2(:,qq)...
        Results_Test_stat_mat_boot_par3(:,qq) Results_Test_stat_mat_boot_par4(:,qq)];
    vect=Test_stat_mat(qq,:)';
    mbar = mean(mat,1);
    vcov = (mat-mbar)'*(mat-mbar)/Nboot;
    stat = vect'*(vcov\vect);
    % show statistic and associated p-value
    disp([stat  1-chi2cdf(stat, 4)])
end
