function interference(design,linkprob)


R = 5000; 
G = 250; N = 25; n = G*N; beta = 1; gamma = .5; %linkprob = .25;

OLS = zeros(R,3); GMM = zeros(R,3); FSTAGE = zeros(R,5); SE_OLS = OLS; SE_GMM = GMM;

parfor r=1:R
    % data generation
    a = randn(N,G); 
    
    if design==0, e =   randn(N,G); end
    if design==1, e = a+randn(N,G); end
    if design==2, e = exp(3*normcdf(a))+randn(N,G); end
    if design==3, e = sin(3*normcdf(a))+randn(N,G); end
    
    x =1+randn(N,G); %>0; 
    
     y = zeros(N,G); Hx = zeros(N,G);
    z1 = zeros(N,G); z2 = zeros(N,G); z3 = zeros(N,G);
    
    A = zeros(N,N,G); 
    H = zeros(N,N,G);
    
    XX = zeros(3,3,G); XY = zeros(3,1,G);
    ZZ = zeros(5,5,G); ZX = zeros(5,3,G); ZY = zeros(5,1,G);
    for g = 1:G
        % network generation
        for i=1:N
            for j=i+1:N
                   A(i,j,g) = a(i,g)+a(j,g)>-sqrt(2)*norminv(linkprob); 
                   A(j,i,g) = A(i,j)           ;
            end
        end  
        H(:,:,g) = diag(1./(max(1,sum(A(:,:,g),2)))) *A(:,:,g); Hx(:,g) = H(:,:,g)*x(:,g);
        % outcome equation
        y(:,g) = x(:,g)*beta + Hx(:,g)*gamma + e(:,g);  
        % instrument creation
        for i=1:N
            Ai = A(:,:,g);      Ai(:,i)=[]; Ai(i,:)=[]; Hi = diag(1./(max(1,sum(Ai,2)))) *Ai; 
            xi = x(:,g); xi(i)  =[];
            
            z1(i,g) = mean(Hi   *xi);
            z2(i,g) = mean(Hi*Hi*xi);
            z3(i,g) = mean(Hi*Hi*Hi*xi); % xi'*(Hi'*(Hi*Hi))*xi/(N-1);
        end 
        
        onesN = ones(N,1);
        % moment matrices for OLS
        XX(:,:,g) = [onesN x(:,g) Hx(:,g)]'*[onesN x(:,g) Hx(:,g)];
        XY(:,:,g) = [onesN x(:,g) Hx(:,g)]'*[y(:,g)              ];
        % moment matrices for GMM
        ZZ(:,:,g) = [onesN x(:,g) z1(:,g) z2(:,g) z3(:,g)]'*[onesN x(:,g) z1(:,g) z2(:,g) z3(:,g)];
        ZX(:,:,g) = [onesN x(:,g) z1(:,g) z2(:,g) z3(:,g)]'*[onesN x(:,g) Hx(:,g)                ];
        ZY(:,:,g) = [onesN x(:,g) z1(:,g) z2(:,g) z3(:,g)]'*[y(:,g)                              ]; 
    end
    
    xx = mean(XX,3)*G/n;
    xy = mean(XY,3)*G/n;
    zx = mean(ZX,3)*G/n;
    zz = mean(ZZ,3)*G/n;
    zy = mean(ZY,3)*G/n;
    
    ols = xx\xy;  gmm = (zx'*(zz\zx))\(zx'*(zz\zy)); fstage = zz\zx;
    
    r_ols = zeros(N,G); omega_ols = zeros(3,3,G);
    r_gmm = zeros(N,G); omega_gmm = zeros(5,5,G);
    
    for g=1:G
        r_ols(:,g) = y(:,g)-[onesN x(:,g) Hx(:,g)]*ols;
        r_gmm(:,g) = y(:,g)-[onesN x(:,g) Hx(:,g)]*gmm;
        
        omega_ols(:,:,g) = [onesN x(:,g) Hx(:,g)                ]'*r_ols(:,g)*r_ols(:,g)'*[onesN x(:,g) Hx(:,g)                ];
        omega_gmm(:,:,g) = [onesN x(:,g) z1(:,g) z2(:,g) z3(:,g)]'*r_gmm(:,g)*r_gmm(:,g)'*[onesN x(:,g) z1(:,g) z2(:,g) z3(:,g)];
    end
    
    O_ols = mean(omega_ols,3)*G/n; asyvar_ols = inv(xx)                        *O_ols*inv(xx)                       ; se_ols = sqrt(diag(asyvar_ols)/n);
    O_gmm = mean(omega_gmm,3)*G/n; asyvar_gmm = inv((zx'*(zz\zx)))*(zx'*inv(zz)*O_gmm*inv(zz)*zx)*inv((zx'*(zz\zx))); se_gmm = sqrt(diag(asyvar_gmm)/n);
    
    OLS(r,:) = ols; SE_OLS(r,:) = se_ols;
    GMM(r,:) = gmm; SE_GMM(r,:) = se_gmm; FSTAGE(r,:) = fstage(:,end);
     
end

[mean(OLS); std(OLS); mean(SE_OLS)],
[mean(GMM); std(GMM); mean(SE_GMM)],

[mean(FSTAGE); std(FSTAGE)],


T = (OLS(:,2:3)-ones(R,1)*[beta gamma])./SE_OLS(:,2:3); 
[mean(T); std(T); mean(abs(T)>=1.96)],

T = (GMM(:,2:3)-ones(R,1)*[beta gamma])./SE_GMM(:,2:3); 
[mean(T); std(T); mean(abs(T)>=1.96)],

