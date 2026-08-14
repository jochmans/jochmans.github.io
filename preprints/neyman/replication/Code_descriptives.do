*** CODE FOR "A Neyman Orthogonalization Approach to the Incidental Parameter Problem", by Bonhomme, Jochmans, Weidner
*** Uses data_research.dta, from Ductor, Fafchamps, Goyal and van der Leij (2014): https://direct.mit.edu/rest/article/96/5/936/58200/Social-Networks-and-Research-Output 

*** This code produces the estimation sample

clear all
set more off

* Enter your path here
cd ""

* use data from Ductor et al
use data_research, clear

* select years
keep if year>=1990 & year<=1999

* net out time effects
bys year: egen meanprod=mean(prod)

su meanprod if year==1999

scalar rmean=r(mean)

replace prod=prod/meanprod*rmean

expand nauthors

gen author=auth1

bys article: replace author=auth2 if _n==2

bys article: replace author=auth3 if _n==3

bys article (author): gen size_team=_N 

* keep team size<=2 
drop if size_team>2

* number of publications
bys author (article): gen nb_pub=_N 

bys author (article): egen nb_pub_own=sum(size_team==1)

gen author_obs=0

bys author (article): replace author_obs=1 if _n==_N

* number of solo publications, by author
tab nb_pub_own if author_obs==1

* keep authors if number of solo publications>=2 - sample selection 
keep if nb_pub_own>=2

* recompute size of the team
bys article (author): replace size_team=_N 

*nb of articles sole authored
count if size_team==1

*nb of articles in teams of size 2 (times 2)
count if size_team==2

*nb of authors
count if author_obs==1

* average journal quality
bys author: egen avprod=mean(prod)
su avprod if author_obs==1, det

* between author variance
su avprod
di r(Var)
su prod
di r(Var)

* journal quality
su prod, det

* number of publications per author
bys author: egen nbprod=sum(1)
su nbprod if author_obs==1, det


 
* prepare data

bys article (author): gen vararticle=(_n==1)

gen sumvararticle=sum(vararticle)

replace article=sumvararticle

bys author (article): gen varauthor=(_n==1)

gen sumvarauthor=sum(varauthor)

replace author=sumvarauthor

keep article author prod size_team

order article author prod size_team

sort article author

outsheet article author prod using data_research, nonames replace 

