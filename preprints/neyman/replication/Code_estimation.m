%%% CODE FOR "A Neyman Orthogonalization Approach to the Incidental Parameter Problem in Likelihood Models", by Bonhomme, Jochmans, Weidner

%%% This code estimates the team production model
%%% It produces the point estimates in Tables 1 and 2 in the paper
%%% It also produces Figures 1 and 2 in the supplement

clear
clc

% Load dataset
load('data_research_processed');

% Loop for cross-fitting

Ncf=100;

Results=zeros(Ncf,15);
Results_0=zeros(Ncf,2);
Results_1=zeros(Ncf,2);
Results_2=zeros(Ncf,2);
Results_3=zeros(Ncf,2);
Results_4=zeros(Ncf,2);
Results_5=zeros(Ncf,2);
Results_6=zeros(Ncf,2);
Results_par0=zeros(Ncf,4);
Results_par1=zeros(Ncf,4);
Results_par2=zeros(Ncf,4);
Results_par3=zeros(Ncf,4);
Results_par4=zeros(Ncf,4);
Results_par5=zeros(Ncf,4);
Results_par6=zeros(Ncf,4);

tic

parfor jcf=1:Ncf

    % seed for reproducibility
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
    % Note: in the code, "alpha" is used for what we call "eta" in the text

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

    % first profite out delta = gamma^(-1)

    delta_loop3=fminbnd(@(a) CES_lik_small(a,full(x2B),full(y2B),alpha_tilde), -500,500);

    Vect=log(y2B)-delta_loop3*log((1/2)*x2B*(alpha_tilde.^(1/delta_loop3)));

    lambda_loop3=mean(Vect);

    sigma2_loop3=mean((Vect-lambda_loop3).^2);

    par_loop3=[lambda_loop3 delta_loop3 log(sigma2_loop3)]';

    % variance in teams of size 1

    var_size1_loop3=mean((log(y1B)-log(x1B*alpha_tilde)).^2);

    par_hat=[1/par_loop3(2) par_loop3(1) exp(par_loop3(3)) var_size1_loop3];

    Results_par0(jcf,:)=par_hat;


    %%% Order q=1 to 6

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

    % minimization of projected score : orders 1 to 6 - to obtain
    % orthogonalized estimators

    warning('off')
    options=optimset('maxiter',5000,'Display','off');

    par_hat1=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,1),par_hat,options);

    % In some samples, the 1st order orthogonal estimate of the 1st parameter may be very close to 0, causing numerical overflow
    % We regularize it using a small number
    if abs(par_hat1(1))<0.01
        par_hat1(1)=sign(par_hat1(1))*0.01;
    end

    par_hat2=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,2),par_hat,options);

    par_hat3=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,3),par_hat,options);

    par_hat4=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,4),par_hat,options);

    par_hat5=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,5),par_hat,options);

    par_hat6=fminunc(@(a) compute_score_f(a,Mat_meanoutput,ones(size(x2B,1),3),Mat_alpha_prelim,6),par_hat,options);

    % Store parameter estimates (for Table 1) - Note: the second parameter is log(beta)

    Results_par1(jcf,:)=par_hat1;
    Results_par2(jcf,:)=par_hat2;
    Results_par3(jcf,:)=par_hat3;
    Results_par4(jcf,:)=par_hat4;
    Results_par5(jcf,:)=par_hat5;
    Results_par6(jcf,:)=par_hat6;

    %%% THIS PART IS FOR FIGURE 2 IN THE SUPPLEMENT
    % Two alternative expressions for beta

    par=par_hat;

    beta_par=exp(par(2));

    beta_model=(mean(exp(par(1)*Mat_meanoutput(:,3)))/((mean(exp(par(1)*Mat_meanoutput(:,1)))+mean(exp(par(1)*Mat_meanoutput(:,2))))/2))^(1/par(1))...
        *exp((par(4)-par(3))*par(1)/2);

    Results_0(jcf,:)=[beta_par beta_model];

    par=par_hat1;

    beta_par=exp(par(2));

    beta_model=(mean(exp(par(1)*Mat_meanoutput(:,3)))/((mean(exp(par(1)*Mat_meanoutput(:,1)))+mean(exp(par(1)*Mat_meanoutput(:,2))))/2))^(1/par(1))...
        *exp((par(4)-par(3))*par(1)/2);

    Results_1(jcf,:)=[beta_par beta_model];

    par=par_hat2;

    beta_par=exp(par(2));

    beta_model=(mean(exp(par(1)*Mat_meanoutput(:,3)))/((mean(exp(par(1)*Mat_meanoutput(:,1)))+mean(exp(par(1)*Mat_meanoutput(:,2))))/2))^(1/par(1))...
        *exp((par(4)-par(3))*par(1)/2);

    Results_2(jcf,:)=[beta_par beta_model];

    par=par_hat3;

    beta_par=exp(par(2));

    beta_model=(mean(exp(par(1)*Mat_meanoutput(:,3)))/((mean(exp(par(1)*Mat_meanoutput(:,1)))+mean(exp(par(1)*Mat_meanoutput(:,2))))/2))^(1/par(1))...
        *exp((par(4)-par(3))*par(1)/2);

    Results_3(jcf,:)=[beta_par beta_model];

    par=par_hat4;

    beta_par=exp(par(2));

    beta_model=(mean(exp(par(1)*Mat_meanoutput(:,3)))/((mean(exp(par(1)*Mat_meanoutput(:,1)))+mean(exp(par(1)*Mat_meanoutput(:,2))))/2))^(1/par(1))...
        *exp((par(4)-par(3))*par(1)/2);

    Results_4(jcf,:)=[beta_par beta_model];

    par=par_hat5;

    beta_par=exp(par(2));

    beta_model=(mean(exp(par(1)*Mat_meanoutput(:,3)))/((mean(exp(par(1)*Mat_meanoutput(:,1)))+mean(exp(par(1)*Mat_meanoutput(:,2))))/2))^(1/par(1))...
        *exp((par(4)-par(3))*par(1)/2);

    Results_5(jcf,:)=[beta_par beta_model];

    par=par_hat6;

    beta_par=exp(par(2));

    beta_model=(mean(exp(par(1)*Mat_meanoutput(:,3)))/((mean(exp(par(1)*Mat_meanoutput(:,1)))+mean(exp(par(1)*Mat_meanoutput(:,2))))/2))^(1/par(1))...
        *exp((par(4)-par(3))*par(1)/2);

    Results_6(jcf,:)=[beta_par beta_model];

    %%% THIS PART IS TO PRODUCE THE ALLOCATION NUMBERS IN TABLE 2
    % mean output, data

    Resinter=zeros(15,1);

    Resinter(1)=mean(exp(Mat_meanoutput(:,3)));

    % mean output, model

    logYsim=par_hat(2)+(par_hat(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(:,1)).^par_hat(1)+exp(Mat_alpha_prelim(:,2)).^par_hat(1))));

    Ysim=exp(logYsim+1/2*par_hat(3));

    Resinter(2)=mean(Ysim);

    % mean output, model, corrected to orders 1-6
    % This is for column 1 in Table 2

    Ysim1full=0;
    Ysim2full=0;
    Ysim3full=0;
    Ysim4full=0;
    Ysim5full=0;
    Ysim6full=0;
    for ii=1:size(x2B,1)
        % correc 1
        [w,SigmaWW,b1,~]=team_CES_order1_Teams_1_2([Mat_meanoutput(ii,1),Mat_meanoutput(ii,2)],[1,1],[],par_hat1(4),[exp(Mat_alpha_prelim(ii,1)),exp(Mat_alpha_prelim(ii,2))],par_hat1(1));
        MM=exp((par_hat1(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii,1)).^par_hat1(1)+exp(Mat_alpha_prelim(ii,2)).^par_hat1(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Ysim1full=Ysim1full+exp(par_hat1(2)+1/2*par_hat1(3))*MM;
        % correc 2
        [w,SigmaWW,b1,~]=team_CES_order2_Teams_1_2([Mat_meanoutput(ii,1),Mat_meanoutput(ii,2)],[1,1],[],par_hat2(4),[exp(Mat_alpha_prelim(ii,1)),exp(Mat_alpha_prelim(ii,2))],par_hat2(1));
        MM=exp((par_hat2(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii,1)).^par_hat2(1)+exp(Mat_alpha_prelim(ii,2)).^par_hat2(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Ysim2full=Ysim2full+exp(par_hat2(2)+1/2*par_hat2(3))*MM;
        % correc 3
        [w,SigmaWW,b1,~]=team_CES_order3_Teams_1_2([Mat_meanoutput(ii,1),Mat_meanoutput(ii,2)],[1,1],[],par_hat3(4),[exp(Mat_alpha_prelim(ii,1)),exp(Mat_alpha_prelim(ii,2))],par_hat3(1));
        MM=exp((par_hat3(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii,1)).^par_hat3(1)+exp(Mat_alpha_prelim(ii,2)).^par_hat3(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Ysim3full=Ysim3full+exp(par_hat3(2)+1/2*par_hat3(3))*MM;
        % correc 4
        [w,SigmaWW,b1,~]=team_CES_order4_Teams_1_2([Mat_meanoutput(ii,1),Mat_meanoutput(ii,2)],[1,1],[],par_hat4(4),[exp(Mat_alpha_prelim(ii,1)),exp(Mat_alpha_prelim(ii,2))],par_hat4(1));
        MM=exp((par_hat4(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii,1)).^par_hat4(1)+exp(Mat_alpha_prelim(ii,2)).^par_hat4(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Ysim4full=Ysim4full+exp(par_hat4(2)+1/2*par_hat4(3))*MM;
        % correc 5
        [w,SigmaWW,b1,~]=team_CES_order5_Teams_1_2([Mat_meanoutput(ii,1),Mat_meanoutput(ii,2)],[1,1],[],par_hat5(4),[exp(Mat_alpha_prelim(ii,1)),exp(Mat_alpha_prelim(ii,2))],par_hat5(1));
        MM=exp((par_hat5(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii,1)).^par_hat5(1)+exp(Mat_alpha_prelim(ii,2)).^par_hat5(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Ysim5full=Ysim5full+exp(par_hat5(2)+1/2*par_hat5(3))*MM;
        % correc 6
        [w,SigmaWW,b1,~]=team_CES_order6_Teams_1_2([Mat_meanoutput(ii,1),Mat_meanoutput(ii,2)],[1,1],[],par_hat6(4),[exp(Mat_alpha_prelim(ii,1)),exp(Mat_alpha_prelim(ii,2))],par_hat6(1));
        MM=exp((par_hat6(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii,1)).^par_hat6(1)+exp(Mat_alpha_prelim(ii,2)).^par_hat6(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Ysim6full=Ysim6full+exp(par_hat6(2)+1/2*par_hat6(3))*MM;
    end

    Ysim1full=Ysim1full/size(x2B,1);
    Ysim2full=Ysim2full/size(x2B,1);
    Ysim3full=Ysim3full/size(x2B,1);
    Ysim4full=Ysim4full/size(x2B,1);
    Ysim5full=Ysim5full/size(x2B,1);
    Ysim6full=Ysim6full/size(x2B,1);

    Resinter(3)=mean(Ysim1full);
    Resinter(4)=mean(Ysim2full);
    Resinter(5)=mean(Ysim3full);
    Resinter(6)=mean(Ysim4full);
    Resinter(7)=mean(Ysim5full);
    Resinter(8)=mean(Ysim6full);


    % mean output in a counterfactual random allocation
    % This is for column 1 in Table 2

    logYrand_inter=par_hat(2)+(par_hat(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(:,1)).^par_hat(1)+(exp(Mat_alpha_prelim(:,2)).^par_hat(1))')));

    Yrand=(sum(sum(exp(logYrand_inter+1/2*par_hat(3)))) - sum(exp(diag(logYrand_inter)+1/2*par_hat(3))))/(size(x2B,1)^2-size(x2B,1));

    Resinter(9)=full(Yrand);

    % random allocation, corrected to orders 1-6

    warning('off')

    Yrand1full=0;
    Yrand2full=0;
    Yrand3full=0;
    Yrand4full=0;
    Yrand5full=0;
    Yrand6full=0;

    % Stochastic approximation to the double sum, using M_ii independent pairs
    M_ii=5000;
    i1_ii=randi(size(x2B,1),M_ii,1);
    i2_ii=randi(size(x2B,1),M_ii,1);
    for mm=1:M_ii
        ii1=i1_ii(mm);
        ii2=i2_ii(mm);
        % correc 1
        [w,SigmaWW,b1,~]=team_CES_order1_Teams_1_2([Mat_meanoutput(ii1,1),Mat_meanoutput(ii2,2)],[1,1],[],par_hat1(4),[exp(Mat_alpha_prelim(ii1,1)),exp(Mat_alpha_prelim(ii2,2))],par_hat1(1));
        MM=exp((par_hat1(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii1,1)).^par_hat1(1)+exp(Mat_alpha_prelim(ii2,2)).^par_hat1(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Yrand1full=Yrand1full+exp(par_hat1(2)+1/2*par_hat1(3))*MM;
        % correc 2
        [w,SigmaWW,b1,~]=team_CES_order2_Teams_1_2([Mat_meanoutput(ii1,1),Mat_meanoutput(ii2,2)],[1,1],[],par_hat2(4),[exp(Mat_alpha_prelim(ii1,1)),exp(Mat_alpha_prelim(ii2,2))],par_hat2(1));
        MM=exp((par_hat2(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii1,1)).^par_hat2(1)+exp(Mat_alpha_prelim(ii2,2)).^par_hat2(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Yrand2full=Yrand2full+exp(par_hat2(2)+1/2*par_hat2(3))*MM;
        % correc 3
        [w,SigmaWW,b1,~]=team_CES_order3_Teams_1_2([Mat_meanoutput(ii1,1),Mat_meanoutput(ii2,2)],[1,1],[],par_hat3(4),[exp(Mat_alpha_prelim(ii1,1)),exp(Mat_alpha_prelim(ii2,2))],par_hat3(1));
        MM=exp((par_hat3(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii1,1)).^par_hat3(1)+exp(Mat_alpha_prelim(ii2,2)).^par_hat3(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Yrand3full=Yrand3full+exp(par_hat3(2)+1/2*par_hat3(3))*MM;
        % correc 4
        [w,SigmaWW,b1,~]=team_CES_order4_Teams_1_2([Mat_meanoutput(ii1,1),Mat_meanoutput(ii2,2)],[1,1],[],par_hat4(4),[exp(Mat_alpha_prelim(ii1,1)),exp(Mat_alpha_prelim(ii2,2))],par_hat4(1));
        MM=exp((par_hat4(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii1,1)).^par_hat4(1)+exp(Mat_alpha_prelim(ii2,2)).^par_hat4(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Yrand4full=Yrand4full+exp(par_hat4(2)+1/2*par_hat4(3))*MM;
        % correc 5
        [w,SigmaWW,b1,~]=team_CES_order5_Teams_1_2([Mat_meanoutput(ii1,1),Mat_meanoutput(ii2,2)],[1,1],[],par_hat5(4),[exp(Mat_alpha_prelim(ii1,1)),exp(Mat_alpha_prelim(ii2,2))],par_hat5(1));
        MM=exp((par_hat5(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii1,1)).^par_hat5(1)+exp(Mat_alpha_prelim(ii2,2)).^par_hat5(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Yrand5full=Yrand5full+exp(par_hat5(2)+1/2*par_hat5(3))*MM;
        % correc 6
        [w,SigmaWW,b1,~]=team_CES_order6_Teams_1_2([Mat_meanoutput(ii1,1),Mat_meanoutput(ii2,2)],[1,1],[],par_hat6(4),[exp(Mat_alpha_prelim(ii1,1)),exp(Mat_alpha_prelim(ii2,2))],par_hat6(1));
        MM=exp((par_hat6(1)^(-1))*log((1/2)*((exp(Mat_alpha_prelim(ii1,1)).^par_hat6(1)+exp(Mat_alpha_prelim(ii2,2)).^par_hat6(1)))));
        MM=MM+w'*(SigmaWW\b1);
        Yrand6full=Yrand6full+exp(par_hat6(2)+1/2*par_hat6(3))*MM;
    end

    Yrand1full=Yrand1full/M_ii;
    Yrand2full=Yrand2full/M_ii;
    Yrand3full=Yrand3full/M_ii;
    Yrand4full=Yrand4full/M_ii;
    Yrand5full=Yrand5full/M_ii;
    Yrand6full=Yrand6full/M_ii;

    Resinter(10)=full(Yrand1full);
    Resinter(11)=full(Yrand2full);
    Resinter(12)=full(Yrand3full);
    Resinter(13)=full(Yrand4full);
    Resinter(14)=full(Yrand5full);
    Resinter(15)=full(Yrand6full);

    % store results on average output (for Table 2)

    Results(jcf,:)=Resinter';

end

toc

%%% TABLE: PARAMETERS, POINT ESTIMATES
%%% Table 1
disp('Table 1')
disp(mean([Results_par0(:,1) exp(Results_par0(:,2)) Results_par0(:,3) Results_par0(:,4)]))
disp(mean([Results_par1(:,1) exp(Results_par1(:,2)) Results_par1(:,3) Results_par1(:,4)]))
disp(mean([Results_par2(:,1) exp(Results_par2(:,2)) Results_par2(:,3) Results_par2(:,4)]))
disp(mean([Results_par3(:,1) exp(Results_par3(:,2)) Results_par3(:,3) Results_par3(:,4)]))
disp(mean([Results_par4(:,1) exp(Results_par4(:,2)) Results_par4(:,3) Results_par4(:,4)]))
disp(mean([Results_par5(:,1) exp(Results_par5(:,2)) Results_par5(:,3) Results_par5(:,4)]))
disp(mean([Results_par6(:,1) exp(Results_par6(:,2)) Results_par6(:,3) Results_par6(:,4)]))

%%% TABLE: AVERAGE OUTPUT, POINT ESTIMATES
%%% Table 2

disp('Table 2')
Results_mean=mean(Results);
disp([Results_mean(2:8)' Results_mean(9:15)'])

disp('mean output in data is')
disp(Results_mean(1))

%%% FIGURE: PRODUCTION FUNCTION
%%% This is Figure 2 in the supplement
eta1=(.1:1:10.1)';
eta2=(.1:1:10.1)';
ga=mean(Results_par6(:,1));
bet=mean(exp(Results_par6(:,2)));
sig2=mean(Results_par6(:,3));

figure

plot(eta2,bet*(1/2*((eta1').^(ga)+eta2.^(ga))).^(1/ga)*exp(1/2*sig2))

axis([0 10 0 30])

xlabel('worker 1')
ylabel('output')

hold off

%%% FIGURE: COMPARING ESTIMATES OF BETA 2
%%% This is Figure 1 in the supplement
markers = {'o', 'd', '^', 's'};

figure

hold on

% Loop through the results and plot each point with different markers
for i = 0:3
    x = mean(eval(['Results_', num2str(i), '(:,1)']));
    y = mean(eval(['Results_', num2str(i), '(:,2)']));

    scatter(x, y, 'Marker', markers{i+1}, 'SizeData', 100); % Set size for visibility

    % Annotate the points with order labels
    if i < 3
        text(x + 0.005, y + 0.005, ['Order ', num2str(i)], ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    else
        text(x + 0.005, y + 0.01, 'Orders 3-6', ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    end
end

% Draw the reference line
plot((1.2:0.001:1.4), (1.2:0.001:1.4), 'k--')

% Set the axis limits
axis([1.2 1.4 1.2 1.4])

% Add x and y labels
xlabel('Parameter Beta')
ylabel('Parameter Beta (Model)')

hold off

% save results
save('Results_est')
