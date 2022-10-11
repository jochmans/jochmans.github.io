function peereffects(design,linkprob)


R = 5000; 
G = 250; N = 25; n = G*N; beta = 1; gamma = .5; delta = .5;

BDF = zeros(R,4); GMM = zeros(R,4); FSTAGE1 = zeros(R,6); FSTAGE2 = FSTAGE1;  SE_BDF = BDF; SE_GMM = GMM;

parfor r=1:R
    % data generation
    a = randn(N,G); 
    
    if design==0, e =   randn(N,G); end
    if design==1, e = a+randn(N,G); end
    if design==2, e = exp(3*normcdf(a))+randn(N,G); end
    if design==3, e = sin(3*normcdf(a))+randn(N,G); end

    
    %if design==4, e = exp(3*betainv(normcdf(a),2,6))+randn(N,G); end
    %if design==5, e = sin(3*betainv(normcdf(a),2,6))+randn(N,G); end
    %if design==6, e = exp(3*betainv(normcdf(a),6,2))+randn(N,G); end
    %if design==7, e = sin(3*betainv(normcdf(a),6,2))+randn(N,G); end
        
    x =1+randn(N,G); %>0;
    %x =randn(N,G)>0;
    
     y = zeros(N,G); Hx = zeros(N,G); Hy = zeros(N,G);
    z1 = zeros(N,G); z2 = zeros(N,G); z3 = zeros(N,G); z4 = zeros(N,G);
    H2x = Hx; H3x = Hx; H4x = Hx;
    
    A = zeros(N,N,G); 
    H = zeros(N,N,G);
    

    ZZb = zeros(6,6,G); ZXb = zeros(6,4,G); ZYb = zeros(6,1,G);
    ZZ  = zeros(6,6,G); ZX  = zeros(6,4,G); ZY  = zeros(6,1,G);
    for g = 1:G
        % network generation
        for i=1:N
            for j=i+1:N
                A(i,j,g) = a(i,g)+a(j,g)>-sqrt(2)*norminv(linkprob); 
                A(j,i,g) = A(i,j);
            end
        end  
        H(:,:,g) = diag(1./(max(1,sum(A(:,:,g),2)))) *A(:,:,g); Hx(:,g) = H(:,:,g)*x(:,g);
        % outcome equation
        y(:,g) = (eye(N)-delta*H(:,:,g))\(x(:,g)*beta + Hx(:,g)*gamma + e(:,g)); Hy(:,g) = H(:,:,g)*y(:,g);
        % instrument creation
        % For BDF
        H2x(:,g) = H(:,:,g)*H(:,:,g)*x(:,g);
        H3x(:,g) = H(:,:,g)*H(:,:,g)*H(:,:,g)*x(:,g);
        H4x(:,g) = H(:,:,g)*H(:,:,g)*H(:,:,g)*H(:,:,g)*x(:,g);
        % For our procedure
        for i=1:N
            Ai = A(:,:,g);      Ai(:,i)=[]; Ai(i,:)=[]; Hi = diag(1./(max(1,sum(Ai,2)))) *Ai; 
            xi = x(:,g); xi(i)  =[];
            
            z1(i,g) = mean(Hi   *xi);
            z2(i,g) = mean(Hi*Hi*xi);
            z3(i,g) = mean(Hi*Hi*Hi*xi); % xi'*(Hi'*Hi)*xi/(N-1);
            z4(i,g) = mean(Hi*Hi*Hi*Hi*xi); % xi'*(Hi'*(Hi*Hi))*xi/(N-1);
        end 
        
        onesN = ones(N,1);
        % moment matrices for BDF
        ZZb(:,:,g) = [onesN x(:,g) Hx(:,g) H2x(:,g) H3x(:,g) H4x(:,g)]'*[onesN x(:,g) Hx(:,g) H2x(:,g) H3x(:,g) H4x(:,g)];
        ZXb(:,:,g) = [onesN x(:,g) Hx(:,g) H2x(:,g) H3x(:,g) H4x(:,g)]'*[onesN x(:,g) Hx(:,g) Hy(:,g)                   ];
        ZYb(:,:,g) = [onesN x(:,g) Hx(:,g) H2x(:,g) H3x(:,g) H4x(:,g)]'*[y(:,g)                                         ]; 
        % moment matrices for GMM
        ZZ(:,:,g) = [onesN x(:,g) z1(:,g) z2(:,g) z3(:,g) z4(:,g)]'*[onesN x(:,g) z1(:,g) z2(:,g) z3(:,g) z4(:,g)];
        ZX(:,:,g) = [onesN x(:,g) z1(:,g) z2(:,g) z3(:,g) z4(:,g)]'*[onesN x(:,g) Hx(:,g) Hy(:,g)                ];
        ZY(:,:,g) = [onesN x(:,g) z1(:,g) z2(:,g) z3(:,g) z4(:,g)]'*[y(:,g)                                      ]; 
    end
    

    zx = mean(ZX,3)*G/n; zxb = mean(ZXb,3)*G/n;
    zz = mean(ZZ,3)*G/n; zzb = mean(ZZb,3)*G/n;
    zy = mean(ZY,3)*G/n; zyb = mean(ZYb,3)*G/n;
    
    bdf = (zxb'*(zzb\zxb))\(zxb'*(zzb\zyb)); 
    gmm = (zx'*(zz\zx))\(zx'*(zz\zy)); fstage = zz\zx;
    
    r_bdf = zeros(N,G); omega_bdf = zeros(6,6,G);
    r_gmm = zeros(N,G); omega_gmm = zeros(6,6,G);
    
    for g=1:G
        r_bdf(:,g) = y(:,g)-[onesN x(:,g) Hx(:,g) Hy(:,g)]*bdf;
        r_gmm(:,g) = y(:,g)-[onesN x(:,g) Hx(:,g) Hy(:,g)]*gmm;
        
        omega_bdf(:,:,g) = [onesN x(:,g) Hx(:,g) H2x(:,g) H3x(:,g) H4x(:,g)]'*r_bdf(:,g)*r_bdf(:,g)'*[onesN x(:,g) Hx(:,g) H2x(:,g) H3x(:,g) H4x(:,g)];
        omega_gmm(:,:,g) = [onesN x(:,g) z1(:,g)  z2(:,g)  z3(:,g)  z4(:,g)]'*r_gmm(:,g)*r_gmm(:,g)'*[onesN x(:,g) z1(:,g)  z2(:,g)  z3(:,g)  z4(:,g)];
    end
    
    O_bdf = mean(omega_bdf,3)*G/n; asyvar_bdf = inv((zxb'*(zzb\zxb)))*(zxb'*inv(zzb)*O_bdf*inv(zzb)*zxb)*inv((zxb'*(zzb\zxb))); se_bdf = sqrt(diag(asyvar_bdf)/n);
    O_gmm = mean(omega_gmm,3)*G/n; asyvar_gmm = inv((zx'*(zz\zx)))*(zx'*inv(zz)*O_gmm*inv(zz)*zx)*inv((zx'*(zz\zx)))          ; se_gmm = sqrt(diag(asyvar_gmm)/n);
    
    BDF(r,:) = bdf; SE_BDF(r,:) = se_bdf;
    GMM(r,:) = gmm; SE_GMM(r,:) = se_gmm; FSTAGE1(r,:) = fstage(:,end-1); FSTAGE2(r,:) = fstage(:,end);
     
end

display('-------------------------------------------------');

design, linkprob,

[mean(BDF); std(BDF); mean(SE_BDF)],
[mean(GMM); std(GMM); mean(SE_GMM)],

[mean(FSTAGE1); std(FSTAGE1)],
[mean(FSTAGE2); std(FSTAGE2)],

T = (BDF(:,2:4)-ones(R,1)*[beta gamma delta])./SE_BDF(:,2:4); 
[mean(T); std(T); mean(abs(T)>=1.96)],

T = (GMM(:,2:4)-ones(R,1)*[beta gamma delta])./SE_GMM(:,2:4); 
[mean(T); std(T); mean(abs(T)>=1.96)],

display('-------------------------------------------------');
display('-------------------------------------------------');

