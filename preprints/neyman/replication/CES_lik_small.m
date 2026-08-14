%%% CODE FOR "A Neyman Orthogonalization Approach to the Incidental Parameter Problem in Likelihood Models", by Bonhomme, Jochmans, Weidner

% Routine to estimate the CES model

function obj=CES_lik_small(par,x,y,alpha_tilde)

Vect=log(y)-par*log((1/2)*x*(alpha_tilde.^(1/par)));

obj=sum((Vect-mean(Vect)).^2);

end
