%%% CODE FOR "A Neyman Orthogonalization Approach to the Incidental Parameter Problem in Likelihood Models", by Bonhomme, Jochmans, Weidner

% This code prepares the dataset for matlab, using the output from "Code_descriptives.do"

clear 
clc

% Load data (generated using Code_descriptives.do)
data=load('data_research.out');

Article=data(:,1);

ID=data(:,2);

Y=data(:,3);

N_res=max(Article);

N_id=max(ID);

[N,~]=size(Article);

% Transform data format

P=sparse(N,N_res);

Ymean=zeros(N_res,1);

for i=1:N_res
    P(:,i)=(Article==i);
    Ymean(i)=mean(Y(Article==i));
    disp(i)
end

Q=sparse(N,N_id);

for j=1:N_id
    Q(:,j)=(ID==j);
    disp(j)
end

% This is an N_res*N_id matrix
A=P'*Q;

% Remove 1 article that is a duplicate
A=min(A,1);

save('data_research_processed')

