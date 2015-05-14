/**************************************/
/* Beginning of first example in text */
/**************************************/

	set more off
	capture log close
	cd "C:\Users\mbaker\Desktop\SJSub-Baker-AMCMC"
	
sjlog using mcmc_p1, replace
clear all
sysuse auto
mata: 
function lregeval(M,todo,b,crit,s,H)
{
real colvector p1, p2
real colvector y1
	p1=moptimize_util_xb(M,b,1)
	p2=moptimize_util_xb(M,b,2)
	y1=moptimize_util_depvar(M,1)
	crit=-(y1:-p1)'*(y1:-p1)/(2*exp(p2))- ///
		rows(y1)/2*p2
}

M=moptimize_init()
moptimize_init_evaluator(M,&lregeval())
moptimize_init_evaluatortype(M,"d0")
moptimize_init_depvar(M,1,"mpg")
moptimize_init_eq_indepvars(M,1,"price weight displacement")
moptimize_init_eq_indepvars(M,2,"")
moptimize(M)
moptimize_result_display(M)
end
sjlog close, replace

/***************************************/
/* Beginning of second example in text */
/***************************************/

sjlog using mcmc_p2, replace
set seed 8675309
mata:
alginfo="moptimize","d0","mwg" 
b_mwg=amcmc(alginfo,&lregeval(),J(1,5,0),
                I(5),10000,50,2/3,.4,
				arate=.,vals=.,lambda=.,.,M)

st_matrix("b_mwg",mean(b_mwg))
st_matrix("V_mwg",variance(b_mwg))
end
local names eq1:price eq1:weight eq1:displacement eq1:_cons eq2:_cons
mat colnames b_mwg=`names'
mat colnames V_mwg=`names'
mat rownames V_mwg=`names'
ereturn post b_mwg V_mwg
ereturn display 
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p3,replace
mata:
arate'
max(vals),mean(vals)
end
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p4,replace
preserve
clear
local varnames price weight displacement constant std_dev
getmata (`varnames')=b_mwg
getmata vals=vals
gen t=_n
local graphs
local tgraphs
foreach var of local varnames {
	quietly {
	histogram `var', saving(`var', replace) nodraw
	twoway line `var' t, saving(t`var', replace) nodraw
	        }
	local graphs "`graphs' `var'.gph"
	local tgraphs "`tgraphs' t`var'.gph"
				}
	histogram vals, saving(vals,replace) nodraw
	twoway line vals t, saving(vals_t,replace) nodraw
				
graph combine `graphs' vals.gph
graph export vals_mwg.eps, replace
graph combine `tgraphs' vals_t.gph
graph export valst_mwg.eps, replace
restore
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p45, replace
mata:
A=amcmc_init()
amcmc_alginfo(A,("global","d0","moptimize"))
amcmc_lnf(A,&lregeval())
amcmc_xinit(A,J(1,5,0))
amcmc_Vinit(A,I(5))
amcmc_model(A,M)
amcmc_draws(A,4000)
amcmc_damper(A,2/3)
amcmc_draw(A)
end
sjlog close, replace

/**************************************/
/* Beginning of third example in text */
/**************************************/

sjlog using mcmc_p5, replace
set seed 8675309
mata: 
alginfo="global","d0","moptimize"
b_glo=amcmc(alginfo,&lregeval(),J(1,5,0),
                I(5),12000,2000,2/3,.4,
				arate=.,vals=.,lambda=.,.,M)

st_matrix("b_glo",mean(b_glo))
st_matrix("V_glo",variance(b_glo))
end
local names eq1:price eq1:weight eq1:displacement eq1:_cons eq2:_cons
mat colnames b_glo=`names'
mat colnames V_glo=`names'
mat rownames V_glo=`names'
ereturn post b_glo V_glo
ereturn display 
sjlog close, replace
preserve
clear
local varnames price weight displacement constant std_dev
getmata (`varnames')=b_glo
getmata vals=vals
gen t=_n
local graphs
local tgraphs
foreach var of local varnames {
	quietly {
	histogram `var', saving(`var', replace) nodraw
	twoway line `var' t, saving(t`var', replace) nodraw
	        }
	local graphs "`graphs' `var'.gph"
	local tgraphs "`tgraphs' t`var'.gph"
				}
	histogram vals, saving(vals,replace) nodraw
	twoway line vals t, saving(vals_t,replace) nodraw
				
graph combine `graphs' vals.gph
graph export vals_glo.eps, replace
graph combine `tgraphs' vals_t.gph
graph export valst_glo.eps, replace
restore
sjlog using mcmc_p6, replace
mata: 
alginfo="mwg","d0","moptimize"
b_start=amcmc(alginfo,&lregeval(),J(1,5,0),
                I(5),5*1000,5*100,2/3,.4,
				arate=.,vals=.,lambda=.,.,M)
alginfo="global","d0","moptimize"
b_glo2=amcmc(alginfo,&lregeval(),mean(b_start),
                variance(b_start),11000,1000,2/3,.4,
				arate=.,vals=.,lambda=.,.,M)
st_matrix("b_glo2",mean(b_glo2))
st_matrix("V_glo2",variance(b_glo2))	
end
			
local names eq1:price eq1:weight eq1:displacement eq1:_cons eq2:_cons
mat colnames b_glo2=`names'
mat colnames V_glo2=`names'
mat rownames V_glo2=`names'
ereturn post b_glo2 V_glo2
ereturn display 
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p65, replace
mata:
alginfo="mwg","d0","moptimize"
b_start=amcmc(alginfo,&lregeval(),J(1,5,0),
                I(5),5*1000,5*100,2/3,.4,
				arate=.,vals=.,lambda=.,.,M)
alginfo="global","d0","moptimize"
b_glo3=amcmc(alginfo,&lregeval(),mean(b_start),
               variance(b_start),10000,0,.,.,
			    arate=.,vals=.,(2.38^2/5),.,M)
arate'
mean(b_glo3)'				
end
sjlog close, replace

/* Code to produce graphs - code not included in paper */

preserve
clear 
getmata (b_glo2*)=b_glo2
getmata vals=vals
gen t=_n
local graphs
local tgraphs
forvalues i=1/5 {
	quietly {
	histogram b_glo2`i', saving(b_glo2`i', replace) nodraw
	twoway line b_glo2`i' t, saving(bt_glo2`i', replace) nodraw
	        }
	local graphs "`graphs' b_glo2`i'.gph"
	local tgraphs "`tgraphs' bt_glo2`i'.gph"
				}
	histogram vals, saving(vals,replace) nodraw
	twoway line vals t, saving(vals_t,replace) nodraw
				
graph combine `graphs' vals.gph
graph export vals_glo2.eps, replace
graph combine `tgraphs' vals_t.gph
graph export valst_glo2.eps, replace

/*************************************/
/* Next example - Censored Quantiles */
/*************************************/

sjlog using mcmc_p66, replace
mata: 
void cqregeval(M,todo,b,crit,g,H) {
	real colvector u,Xb,y,C
	real scalar tau
	
	Xb    =moptimize_util_xb(M,b,1)		 
	y     =moptimize_util_depvar(M,1)	
	tau   =moptimize_util_userinfo(M,1) 	
	C     =moptimize_util_userinfo(M,2)	
	u     =(y:-rowmax((C,Xb)))
	crit  =-colsum(u:*(tau:-(u:<0))) 
}
end
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p67, replace
webuse laborsub, clear
gen censorpoint=0
mata: 
M=moptimize_init()
moptimize_init_evaluator(M,&cqregeval())
moptimize_init_depvar(M,1,"whrs")
moptimize_init_eq_indepvars(M,1,"kl6 k618 wa")
tau=.6
moptimize_init_userinfo(M,1,tau)
st_view(C=.,.,"censorpoint")
moptimize_init_userinfo(M,2,C)
moptimize_init_evaluatortype(M,"d0")
end
sjlog close, replace

mata: moptimize_evaluate(M)

/* Code snippet */

sjlog using mcmc_p68, replace
mata: 
alginfo="mwg","d0","moptimize"
b_start=amcmc(alginfo,&cqregeval(),J(1,4,0),
                I(4),5000,1000,2/3,.4,
				arate=.,vals=.,lambda=.,.,M)
alginfo="global","d0","moptimize"
b_end=amcmc(alginfo,&cqregeval(),mean(b_start),
               variance(b_start),20000,10000,1,.234,
			    arate=.,vals=.,lambda=.,.,M)
end
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p69, replace
set seed 584937
qui mcmccqreg whrs kl6 k618 wa, tau(.6) sampler("mwg") draws(5000) ///
	burn(1000) dampparm(.667) arate(.4) censorvar(censorpoint)
	mat binit=e(b)
	mat V=e(V)
mcmccqreg whrs kl6 k618 wa, tau(.6) sampler("global") draws(20000) ///
	burn(10000) arate(.234) saving(lsub_draws) replace ///
	from(binit) fromv(V)
sjlog close, replace

/***************************************/
/* Beginning of fourth example in text */
/***************************************/

clear all
sjlog using mcmc_p7, replace
set seed 262728
mata:
real scalar ln_fun(x)
{
	return(-x[1]^2-1/2*x[2]^2+x[1]*x[2]-.05*(x[3]-100)^2)
}
B=(1,1,0) \ (0,0,1)
alginfo="standalone","block"
x_block=amcmc(alginfo,&ln_fun(),J(1,3,0),
                I(3),4000,200,2/3,.4,
				arate=.,vals=.,lambda=.,B)
end
sjlog close, replace

/* Additional code to make graphs - not in text */

local graphs
local tgraphs

getmata (x_*)=x_block
getmata vals
gen t=_n
forvalues i=1/3 {
	quietly {
	histogram x_`i', saving(x_block`i', replace) nodraw
	twoway line x_`i' t, saving(xt_block`i', replace) nodraw
	        }
	local graphs "`graphs' x_block`i'.gph"
	local tgraphs "`tgraphs' xt_block`i'.gph"
				}
	histogram vals, saving(vals,replace) nodraw
	twoway line vals t, saving(vals_t,replace) nodraw
				
graph combine `graphs' vals.gph
graph export x_block.eps, replace
graph combine `tgraphs' vals_t.gph
graph export xt_block.eps, replace

/* Code snippet */

sjlog using mcmc_p75, replace
mata:
A=amcmc_init()
amcmc_lnf(A,&ln_fun())
amcmc_alginfo(A,("standalone","block"))
amcmc_draws(A,4000)
amcmc_burn(A,200)
amcmc_damper(A,2/3)
amcmc_xinit(A,J(1,3,0))
amcmc_Vinit(A,I(3))
amcmc_blocks(A,B)
amcmc_draw(A)
end
sjlog close, replace

/********************************************/
/* Discussion of Bayesian Mixed Logit Model */
/********************************************/

cd "C:\Users\mbaker\Desktop\SJSub-Baker-AMCMC"

/* Code snippet */

sjlog using mcmc_p9, replace
clear all
set more off
use http://fmwww.bc.edu/repec/bocode/t/traindata.dta

set seed 90210
mixlogit y, rand(price contract local wknown tod seasonal) group(gid) id(pid)
sjlog close, replace

/* Next snippet */

sjlog using mcmc_p91, replace
mata:
real matrix drawb_betaW(beta,W) {
	return(mean(beta)+rnormal(1,cols(beta),0,1)*cholesky(W)')
}
end
sjlog close, replace

/* Next snippet */

sjlog using mcmc_p92, replace
mata
real matrix drawW_bbeta(beta,b) 
{
    v=rnormal(cols(b)+rows(beta),cols(b),0,1)
	S1=variance(beta)
	S=invsym((cols(b)*I(cols(b))+rows(beta)*S1)/(cols(b)+rows(beta)))
	L=cholesky(S)
	R=(L*v')*(L*v')'/(cols(b)+rows(beta))
	return(invsym(R))
}	
end
sjlog close, replace

/* Next snippet */

sjlog using mcmc_p93, replace
mata:
st_view(y=.,.,"y")
st_view(X=.,.,"price contract local wknown tod seasonal")
st_view(pid=.,.,"pid")
st_view(gid=.,.,"gid")
end
sjlog close, replace

/* Next snippet */

sjlog using mcmc_p94, replace
mata:
real scalar lnbetan_bW(betaj,b,W,yj,Xj)
{
	Uj=rowsum(Xj:*betaj)
	Uj=colshape(Uj,4)
	lnpj=rowsum(Uj:*colshape(yj,4)):-
		ln(rowsum(exp(Uj)))
	var=-1/2*(betaj:-b)*invsym(W)*(betaj:-b)'-
		1/2*ln(det(W))-cols(betaj)/2*ln(2*pi())
	llj=var+sum(lnpj)
	return(llj)
}
end
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p95, replace
mata
m=panelsetup(pid,1)
Ap=amcmc_init()
amcmc_damper(Ap,1)
amcmc_alginfo(Ap,("standalone","global"))
amcmc_append(Ap,"overwrite")
amcmc_lnf(Ap,&lnbetan_bW())
amcmc_draws(Ap,1)
amcmc_append(Ap,"overwrite")
amcmc_reeval(Ap,"reeval")
A=J(rows(m),1,Ap)
end
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p96, replace
mata
Args=J(rows(m),4,NULL)
b=J(1,6,0)
W=I(6)*6
beta=b:+sqrt(diagonal(W))':*rnormal(rows(m),cols(b),0,1)
for (i=1;i<=rows(m);i++) {
	Args[i,1]=&b
	Args[i,2]=&W
	Args[i,3]=&panelsubmatrix(y,i,m)
	Args[i,4]=&panelsubmatrix(X,i,m)
	amcmc_args(A[i],Args[i,])
	amcmc_xinit(A[i],b)
	amcmc_Vinit(A[i],W)
					     }
end
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p97, replace
mata
its=20000
burn=10000
bvals=J(0,cols(beta),.)
Wvals=J(0,cols(rowshape(W,1)),.)
for (i=1;i<=its;i++) {
	b=drawb_betaW(beta,W/rows(m))
	W=drawW_bbeta(beta,b)
	bvals=bvals\b
	Wvals=Wvals\rowshape(W,1)
	beta_old=beta
	for (j=1;j<=rows(A);j++) {
		amcmc_draw(A[j])
		beta[j,]=amcmc_results_lastdraw(A[j])
	}
}
		 
end
sjlog close, replace

/* Code snippet */
/* Unused in the paper!
sjlog using mcmc_p98, replace
mata:
bvals=bvals[burn::its,]
Wvals=Wvals[burn::its,]
Kept=colshape(1::10000,10)
bKept=bvals[Kept[,10],]
WKept=Wvals[Kept[,10],]
mean(bKept)', sqrt(diagonal(variance(bKept)))
end
sjlog close, replace

/* Code snippet */

sjlog using mcmc_p99, replace
mata: round(cholesky(colshape(mean(WKept),6)),1e-3)
sjlog close, replace
*/
/* Code snippet */

sjlog using mcmc_p100, replace
set seed 475446
bayesmixedlogit y, rand(price contract local wknown tod seasonal) ///
	group(gid) id(pid) draws(20000) burn(10000) ///
	samplerrand("global") saving(train_draws) replace
sjlog close, replace

/***********************************************/
/* Code snippet describing the basic algorithm */
/***********************************************/

clear*

sjlog using mcmc_p8, replace
mata:
real matrix amcmc_global(f,xinit,Vinit,draws,burn,damper,
					   aopt,arate,val,lam)
{
	real scalar nb,old,pro,i,alpha
	real rowvector xold,xpro,mu
	real matrix Accept,accept,xs,V,Vsq,Vold

	nb=cols(xinit)  /* Initialization */
	xold=xinit      
	lam=2.38^2/nb   
	old=(*f)(xold)     
	val=old
	
	Accept=0            
	xs=xold             
	mu=xold             
	V=Vinit             
	Vold=I(cols(xold))  

	for (i=1;i<=draws;i++)  {
		accept=0		    
		Vsq=cholesky(V)'    /* Prep V for drawing */
		if (hasmissing(Vsq)) {    
			Vsq=cholesky(Vold)'  
			V=Vold
        }
							 
		xpro=xold+lam*rnormal(1,nb,0,1)*Vsq  /* Draw, value calc. */
											
		pro=(*f)(xpro)   
		
		if      (pro==. ) alpha=0    	/* calc. of accept. prob */      
		else if (pro>old) alpha=1    
		else alpha=exp(pro-old)
			
		if (runiform(1,1)<alpha) {
			old=pro
			xold=xpro
			accept=1
        }
									 
		lam=lam*exp(1/(i+1)^damper*(alpha-aopt)) /*update*/
		xs=xs\xold
		val=val\old
		Accept=Accept\accept
		mu=mu+1/(i+1)^damper*(xold-mu)
		Vold=V
		V=V+1/(i+1)^damper*((xold-mu)'(xold-mu)-V)
		_makesymmetric(V)
	}	
					   
	val  =val[burn+1::draws,]
	arate=mean(Accept[burn+1::draws,])
	return(xs[burn+1::draws,])
}
end
sjlog close, replace










