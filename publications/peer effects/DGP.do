clear
clear mata
program drop _all
set obs 5000


forvalues i=1(1)3 {
	gen x`i'=rnormal()+1
}


gen G0=(round(uniform()*200)+1)
sort G0
qui sum G0
local max=r(max)

mata: ivnet_A=asarray_create("real",1)
mata: a=asarray_create("real",1)
mata: x=asarray_create("real",1)
mata: yfin=J(0,1,.)
mata: G=`max'

forvalues i=1(1)`max' {
	qui putmata x0=(x*) if G==`i', replace
	mata: asarray(x,`i',x0)
	mata: U=J(rows(x0),rows(x0),0)
	mata: asarray(ivnet_A,`i',U)
	mata: asarray(a,`i',rnormal(rows(U),1,0,1))
}

mata {

	/// Create Matrix

	Ng=0; beta = 1; gamma = .5; delta = .5;
	for (g=1;g<=G;g++) {
		///Network generation
		Ag=asarray(ivnet_A,g)
		N=rows(Ag)
		Ng=Ng+N
		ag=asarray(a,g)
		for (i=1;i<=N;i++) {
			for (j=i+1;j<=N;j++) {
				Ag[i,j]=(ag[i,]+ag[j,])>-sqrt(2)*invnormal(0.25)
				Ag[j,i]=Ag[i,j]
			}

		}
		asarray(ivnet_A,g,Ag)
	}


	for (g=1;g<=G;g++) {
	Ag=asarray(ivnet_A,g)
	ag=asarray(a,g)
	yg=J(rows(Ag),0,.)
	eg=rnormal(rows(Ag),1,0,1)
	///:+exp(3*normal(ag))
	Hxreg=0
	
		for (k=1;k<=cols(x0);k++) {
			xg=asarray(x,g)[,k]
			Hg = Ag:/rowmax((J(rows(Ag),1,1),rowsum(Ag))); Hx=Hg*xg;
			Hxreg=Hxreg:+Hx*gamma
		}

		yg=lusolve((I(rows(Ag))-delta*Hg),(asarray(x,g)*J(cols(x0),1,1):+Hxreg+eg)); Hy=Hg*yg;
		yfin=(yfin\yg)

	}
	}

getmata y=yfin	
mata: mata matsave ivnet_A ivnet_A, replace 



	