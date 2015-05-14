/* 5-30-2013 - a fix of a bug to the evaluation function */
/* for the moptimize case was repaired.                  */
/* additional helpfile added under the name of amcmc     */
mata: 
mata clear
mata set matastrict on
real matrix amcmc(transmorphic alginfo,    
					  pointer (real scalar function) scalar f,
	                  real rowvector binit,
					  real matrix Vinit,
					  real scalar draws,
					  real scalar burn,
					  real scalar damper,
					  real scalar aopt,
					  transmorphic arate,
					  transmorphic val,
					  transmorphic lam,
					  real matrix blocks,
					  | transmorphic M,
					  string scalar noisy)
{
	string matrix setinf
	string scalar type,opttype,etype,fast,inine
	
	  /* Parse options */
	
	  setinf=amcmc_parser(alginfo)
	  type   =setinf[1]
	  opttype=setinf[2]
	  etype  =setinf[3]
	  fast   =setinf[4]
	  inine  =setinf[5]

	  /* Parse function options */
	
	  if (args()==13 & M=="noisy") {
			M=.
			noisy="noisy"
								   }

	  if (type=="mwg") return(amcmc_mwg(f,binit,Vinit,draws,burn,
	                                    damper,aopt,arate,val,lam,
										M,opttype,etype,inine,noisy))
					  
	  else if (type=="block") return(amcmc_block(f,binit,Vinit,blocks,
	                                             draws,burn,damper,aopt,
					                             arate,val,lam,M,opttype,
					                             etype,fast,inine,noisy))
					  
	  else return(amcmc_global(f,binit,Vinit,draws,burn,damper,aopt,
					           arate,val,lam,M,opttype,etype,fast,
							   inine,noisy))
} 
real matrix amcmc_mwg(pointer (real scalar function) scalar f,
	                  real rowvector xinit,
					  real matrix Vinit,
					  real scalar draws,
					  real scalar burn,
					  real scalar damper,
					  real scalar aopt,
					  transmorphic arate,
					  transmorphic val,
					  transmorphic lam,
					  transmorphic M,
					  string scalar opttype,
					  string scalar etype,
					  string scalar inine,
					  string scalar noisy)
{
	real scalar nb,old,pro,i,j,alpha,its
	real rowvector xold,xpro,sigma
	real matrix s,H,Accept,accept,xs

	/* Initial information */ 
	
	nb   =cols(xinit)
	xold =xinit
	sigma=diag(diagonal(Vinit)) /* Fixed throughout - lambda adjusts */
	
	if (damper!=.) lam=J(1,nb,2.38^2/nb)	/* Set adj. if damper!=. */
	
	old=amcmc_feval(f,xold,M,etype,opttype) 
	val=old			
	
	Accept=J(1,nb,0)
	alpha =J(1,nb,0)
	xs    =xold
	its   =1
	
	for (i=1;i<=draws;i++)  {
		accept=J(1,nb,0)
		for (j=1;j<=nb;j++) {
			xpro=xold
			xpro[j]=xold[j]+rnormal(1,1,0,1)*sqrt(sigma[j,j])*lam[j]
			pro=amcmc_feval(f,xpro,M,etype,opttype)

			if      (pro==.)  alpha[j]=0	// Compute prop. accept. prob.
			else if (pro>old) alpha[j]=1
			else alpha[j]=exp(pro-old)

			if (runiform(1,1)<alpha[j]) {	// Accept proposal or not
				old=pro
				xold=xpro
				accept[j]=1
									}
								
			amcmc_makenoise(its,noisy,old)
			its=its+1
		                }	// Termination of j-loop
						
	    if (damper!=.) lam=lam:*exp(1/i^(damper)*(alpha:-aopt))
		xs  =    xs\xold
		val   =   val\old
		Accept=Accept\accept
					   }
					   
	val  =val[burn+1::draws,]		  // Reporting
	arate=mean(Accept[burn+1::draws,])
	return(xs[burn+1::draws,])
}
real matrix amcmc_global(pointer (real scalar function) scalar f,
	                     real rowvector xinit,
					     real matrix Vinit,
					     real scalar draws,
					     real scalar burn,
					     real scalar damper,
					     real scalar aopt,
					     transmorphic arate,
					     transmorphic val,
					     transmorphic lam,
					     transmorphic M,
					     string scalar opttype,
					     string scalar etype,
					     string scalar fast,
					     string scalar inine,
					     string scalar noisy)
{
	real scalar nb,old,pro,i,alpha
	real rowvector xold,xpro,mu
	real matrix s,H,Accept,accept,xs,V,Vsq,Vold
	
	nb  =cols(xinit)
	xold=xinit

	if (damper!=.) lam=2.38^2/nb	// Set adj. parm. if damper given
	
	old=amcmc_feval(f,xold,M,etype,opttype)	
	val=old
	
	Accept=0
	xs    =xold
	mu    =xold
	V     =Vinit
	Vold  =I(cols(xold))		// Fall back V if V doesn't work	
	V=amcmc_Vsetup(V,damper,fast)
	
	for (i=1;i<=draws;i++)  {
		accept=0
		xpro=xold
		Vsq=amcmc_Vdrawprep(V,Vold,damper,fast)

		if (fast=="fast") xpro=xold:+lam*rnormal(1,nb,0,1):*Vsq
		else              xpro=xold:+lam*rnormal(1,nb,0,1)*Vsq

		pro=amcmc_feval(f,xpro,M,etype,opttype)
	
		if      (pro==.)  alpha=0
		else if (pro>old) alpha=1
		else              alpha=exp(pro-old)

		if (runiform(1,1)<alpha) {
				old=pro
				xold=xpro
				accept=1
								  }
		
		xs=xs\xold
		val=val\old
		Accept=Accept\accept

		amcmc_drawinf_update(mu,xold,V,Vold,lam,alpha,aopt,i,damper,fast)
		amcmc_makenoise(i,noisy,old)
					   }	

	val  =val[burn+1::draws,]
	arate=mean(Accept[burn+1::draws,])
	return(xs[burn+1::draws,])
}
real matrix amcmc_block(pointer (real scalar function) scalar f,
	                    real rowvector xinit,
					    real matrix Vinit,
					    real matrix block,
					    real scalar draws,
					    real scalar burn,
					    real scalar damper,
					    real scalar aopt,
					    transmorphic arate,
					    transmorphic val,
					    transmorphic lam,
					    transmorphic M,
					    string scalar opttype,
					    string scalar etype,
					    string scalar fast,
					    string scalar inine,
					    string scalar noisy)
{
	real scalar nb,old,pro,i,j,alpha,its
	real rowvector xold,xpro,mu
	real matrix s,H,Accept,accept,xs,V,Vsq,Vold

	nb  =rows(block)
	xold=xinit
	its=1
	if (damper!=.) lam=J(1,rows(block),1)*2.38^2/nb

	old=amcmc_feval(f,xold,M,etype,opttype)
	val=old
	
	Accept=J(1,nb,0)
	accept=J(1,nb,0)
	alpha =J(1,rows(block),1)
	xs    =xold
	mu    =xold
	V     =Vinit
	Vold=I(cols(xold))
	V=amcmc_Vsetup(V,damper,fast)

	for (i=1;i<=draws;i++)  {
		accept=J(1,nb,0)

	for (j=1;j<=rows(block);j++) {
		xpro=xold
		Vsq=amcmc_Vdrawprep(V,Vold,damper,fast)

		if (fast=="fast") xpro=xold+block[j,]:*(rnormal(1,cols(xpro),0,1):*Vsq*lam[j])
		else              xpro=xold+block[j,]:*(rnormal(1,cols(xpro),0,1)*Vsq*lam[j])
		
						  
		pro=amcmc_feval(f,xpro,M,etype,opttype)	
		
		if      (pro==.)  alpha[j]=0
		else if (pro>old) alpha[j]=1
		else              alpha[j]=exp(pro-old)
		
		if (runiform(1,1)<alpha[j]) {
				old=pro
				xold=xpro
				accept[j]=1
									}

		amcmc_makenoise(its,noisy,old)
		its=its+1									
									} // Close of j-loop	
		xs   =xs\xold
		val   =val\old
		Accept=Accept\accept
		
		amcmc_drawinf_update(mu,xold,V,Vold,lam,alpha,aopt,i,damper,fast)
			}
					   
	val  =val[burn+1::draws,]
	arate=mean(Accept[burn+1::draws,])
	return(xs[burn+1::draws,])
}
/* Structured AMCMC routines and definitions */
struct amcmc_struct {
	  pointer (real scalar function) scalar f
      real rowvector xinit
	  real matrix Vinit
	  real scalar draws,burn,damper,aopt,passes,totaldraws
	  real matrix val
	  real matrix blocks
	  real matrix lam
	  transmorphic M
	  string scalar noisy
	  string scalar append
	  string scalar reeval
	  string matrix alginfo
	  pointer matrix args
	  
	  /* Returned information */
	  real matrix xvals
	  real matrix mu,Vval
	  real matrix vals
	  real matrix Accept
	  real matrix arate
}
void amcmc_draw(struct amcmc_struct A) 
{
	string matrix alginfo
	// Rearrange and overwrite the user-supplied alg. information
	alginfo=A.alginfo
	A.alginfo=amcmc_parser(alginfo)
	
	if      (A.alginfo[1]=="mwg")   amcmc_inline_mwg(A)
	else if (A.alginfo[1]=="block") amcmc_inline_block(A)
	else                            amcmc_inline_global(A)
}
/* Drawing routines */
void amcmc_inline_global(struct amcmc_struct A)
{
	real scalar draws,aopt,burn,damper,totaldraws,
	            i,old,pro,nb,accept,alpha,itouse,arate
	real matrix xold,xs,V,Vold,val,lam,
	            mu,Accept,xpro,Vsq
	string matrix alginfo
	string scalar etype,opttype,fast,noisy
    pointer (real scalar function) scalar f	
    transmorphic M
    
	draws     =A.draws
	aopt      =A.aopt
	burn      =A.burn
	damper    =A.damper	
	f         =A.f
	M         =A.M
	totaldraws=A.totaldraws
	alginfo   =A.alginfo
	noisy     =A.noisy

	  opttype =alginfo[2]
	  etype   =alginfo[3]
	  fast    =alginfo[4]

	nb=cols(A.xinit) 

    if (A.passes==0) {		/* If the first pass, initialize everything */
	   xold=A.xinit			// User supplied initial xinit
	   V   =A.Vinit			// User supplied initial Vinit
	   Vold  =I(cols(xold))	// backup V 	   
	   if (damper!=.) lam=2.38^2/nb	// Set adj. parm. if damper given
	     else lam=A.lam       // lambda given
	   old=amcmc_feval(f,xold,M,etype,opttype) // evaluation
       val=old              // initial value
	   xs    =J(0,cols(xold),0) // collect the xs.
	   mu    =A.mu          // Set by initialization to xinit
	   A.vals=val           // Hold initial value
	   itouse=0             // No continuing iteration
	   Accept=J(0,1,0)
			}
	else {			       /* If not the first pass, pull in last values */
		xold =A.xvals[rows(A.xvals),.] 
		Vold =A.Vval
		V=Vold		
		mu   =A.mu
		if (A.reeval=="reeval") old=amcmc_feval(f,xold,M,etype,opttype)
		else old  =A.vals[rows(A.vals),.]
		itouse=totaldraws
		lam  =A.lam		
		arate=A.arate	// All the previous are as-is if append or not   
		if (A.append=="append") {
			xs    =A.xvals
		    val   =A.vals			
	        Accept=A.Accept		
			}
		else                    {
		    xs=J(0,cols(xold),0)	
			val=old
			Accept=J(0,1,0)
			}
		  }

	V=amcmc_Vsetup(V,damper,fast)
					   
	for (i=1+itouse;i<=itouse+draws;i++)  {
		accept=0
		xpro=xold

		Vsq=amcmc_Vdrawprep(V,Vold,damper,fast)

		if (fast=="fast") xpro=xold:+lam*rnormal(1,nb,0,1):*Vsq
		else              xpro=xold:+lam*rnormal(1,nb,0,1)*Vsq

	/* Calculate the new objective function value */

		pro=amcmc_feval(f,xpro,M,etype,opttype)
		
		if      (pro==.)  alpha=0
		else if (pro>old) alpha=1
		else              alpha=exp(pro-old)

		if (runiform(1,1)<alpha) {
				old=pro
				xold=xpro
				accept=1
								  }

		xs=xs\xold
		val=val\old
		Accept=Accept\accept
		amcmc_drawinf_update(mu,xold,V,Vold,lam,alpha,aopt,i,damper,fast)
		amcmc_makenoise(i,noisy,old)
						}	

		A.Accept=Accept
		A.lam=lam
		A.Vval=V
		A.mu=mu
	if (A.passes==0) {
		A.xvals=xs[burn+1::rows(xs),]
		A.vals=val[burn+1::rows(xs),]
					}
	else {
		 A.xvals=xs
		 A.vals=val
		 }

	if (A.append=="append" | A.passes==0) A.arate=mean(Accept)
	else A.arate=arate*totaldraws/(totaldraws+draws)+
	             draws/(totaldraws+draws)*mean(Accept)
		 
		A.passes=A.passes+1
		A.totaldraws=A.totaldraws+draws
}
void amcmc_inline_mwg(struct amcmc_struct A)
{
	real scalar draws,aopt,burn,damper,totaldraws,nb,old,pro,i,j,alpha,its,itouse
	real matrix xold,xs,Vold,val,lam,Accept,xpro,sigma,accept,arate
	string matrix alginfo
	string scalar etype,opttype,noisy
	pointer (real scalar function) scalar f
	transmorphic M

	draws=A.draws
	aopt=A.aopt
	burn=A.burn
	damper=A.damper
	f=A.f
	M=A.M
	totaldraws=A.totaldraws
	alginfo=A.alginfo
	noisy=A.noisy
	
		opttype=alginfo[2]
		etype=alginfo[3]
	
	nb=cols(A.xinit)

	if (A.passes==0) {
		xold=A.xinit
		Vold=A.Vinit
		if (damper!=.) lam=J(1,nb,2.38^2/nb)
		else lam=A.lam
		old=amcmc_feval(f,xold,M,etype,opttype)
		val=old
		xs=J(0,nb,0)
		Vold=I(nb)
		A.vals=val
		sigma=diag(diagonal(A.Vinit))
		itouse=0
		Accept=J(0,nb,0)
		             }
	else {
	    xold =A.xvals[rows(A.xvals),.]
		Vold =A.Vval
		sigma=diag(diagonal(Vold))
		if (A.reeval=="reeval") old=amcmc_feval(f,xold,M,etype,opttype)
		else old  =A.vals[rows(A.vals),.]
		itouse=totaldraws
		lam  =A.lam
		arate=A.arate 
		if (A.append=="append") {
			xs=A.xvals
			val=A.vals
			Accept=A.Accept
			}
		else {
			xs=J(0,cols(xold),0)
			val=old
			Accept=J(0,nb,0)
			}
		}

	alpha =J(1,nb,0)
	its   =1

	for (i=1+itouse;i<=itouse+draws;i++)  {
		accept=J(1,nb,0)
	for (j=1;j<=nb;j++) {
		xpro=xold
		xpro[j]=xold[j]+rnormal(1,1,0,1)*sqrt(sigma[j,j])*lam[j]
		
		pro=amcmc_feval(f,xpro,M,etype,opttype)
		
		if      (pro==.)  alpha[j]=0	// Compute prop. accept. prob.
		else if (pro>old) alpha[j]=1
		else alpha[j]=exp(pro-old)

		if (runiform(1,1)<alpha[j]) {	// Accept proposal or not
				old=pro
				xold=xpro
				accept[j]=1
									}
									
		amcmc_makenoise(its,noisy,old)
		its=its+1
		                }	// Termination of j-loop
						
	    if (damper!=.) lam=lam:*exp(1/i^(damper)*(alpha:-aopt))
		  xs  =xs\xold
		  val =val\old
		  Accept=Accept\accept
					   }
				   
	A.Accept=Accept
	A.lam=lam
	A.Vval=sigma
	if (A.passes==0) {
		A.xvals=xs[burn+1::rows(xs),]
		A.vals =val[burn+1::rows(xs),]
		             }
	else {
		A.xvals=xs
		A.vals=val
		 }
	if (A.append=="append" | A.passes==0) A.arate=mean(Accept)
	else A.arate=arate*totaldraws/(totaldraws+draws)+
	       draws/(totaldraws+draws)*mean(Accept)
		 
	A.passes=A.passes+1
	A.totaldraws=A.totaldraws+its
}
void amcmc_inline_block(struct amcmc_struct A)
{
	real scalar draws,aopt,burn,damper,totaldraws,nb,old,pro,i,j,alpha,its,arate,itouse
	real matrix xold,xs,V,Vold,Vsq,val,lam,mu,Accept,xpro,accept,block
	string matrix alginfo
	string scalar etype,opttype,fast,noisy
	pointer (real scalar function) scalar f
	transmorphic M

	draws=A.draws
	aopt=A.aopt
	burn=A.burn
	damper=A.damper
	f=A.f
	M=A.M
	totaldraws=A.totaldraws
	alginfo=A.alginfo
	noisy=A.noisy
	
		opttype=alginfo[2]
		etype=alginfo[3]
		fast=alginfo[4]
	
	nb=rows(A.blocks)
	block=A.blocks

	if (A.passes==0) {
		xold=A.xinit
		V=A.Vinit
		Vold=I(cols(xold))
		if (damper!=.) lam=J(1,nb,2.38^2/nb)
		else lam=A.lam
		old=amcmc_feval(f,xold,M,etype,opttype)
		val=old
		xs=J(0,cols(xold),0)
		V=Vold
		mu=A.mu
		A.vals=val
		Accept=J(0,nb,0)
		itouse=0
		             }
	else {
	    xold =A.xvals[rows(A.xvals),.]
		Vold =A.Vval
		V    =Vold
		mu   =A.mu
		if (A.reeval=="reeval") old=amcmc_feval(f,xold,M,etype,opttype)
		else old  =A.vals[rows(A.vals),.]
	    itouse=totaldraws
		lam  =A.lam
		arate=A.arate
		if (A.append=="append") {
			xs=A.xvals
			val=A.vals
			Accept=A.Accept
	        }
		else {
		    xs=J(0,cols(xold),0)
			val=old
			Accept=J(0,nb,0)
			 }
		}

	V=amcmc_Vsetup(V,damper,fast)

	alpha =J(1,nb,0)
	its   =1
	for (i=1+itouse;i<=itouse+draws;i++)  {
		accept=J(1,nb,0)
	for (j=1;j<=nb;j++) {
		xpro=xold
		Vsq=amcmc_Vdrawprep(V,Vold,damper,fast)

		if (fast=="fast") xpro=xold+block[j,]:*(rnormal(1,cols(xpro),0,1):*Vsq*lam[j])
		else              xpro=xold+block[j,]:*(rnormal(1,cols(xpro),0,1)*Vsq*lam[j])


		pro=amcmc_feval(f,xpro,M,etype,opttype)
		
		if      (pro==.)  alpha[j]=0	// Compute prop. accept. prob.
		else if (pro>old) alpha[j]=1
		else alpha[j]=exp(pro-old)

		if (runiform(1,1)<alpha[j]) {	// Accept proposal or not
				old=pro
				xold=xpro
				accept[j]=1
									}

		amcmc_makenoise(its,noisy,old)
			its=its+1
		                }	// Termination of j-loop

 		  amcmc_drawinf_update(mu,xold,V,Vold,lam,alpha,aopt,i,damper,fast)
		  xs  =xs\xold
		  val =val\old
		  Accept=Accept\accept
					   }

	A.Accept=Accept
	A.lam=lam
	A.Vval=V
	A.mu=mu
	if (A.passes==0) {
		A.xvals=xs[burn+1::rows(xs),]
		A.vals =val[burn+1::rows(xs),]
		             }
	else {
		A.xvals=xs
		A.vals=val
		 }
    if (A.append=="append" | A.passes==0) A.arate=mean(Accept)
	else A.arate=arate*totaldraws/(totaldraws+draws):+
	   draws/(totaldraws+draws)*mean(Accept)
	   
	A.passes=A.passes+1
	A.totaldraws=A.totaldraws+its
}
				
struct amcmc_struct amcmc_init()
{
	struct amcmc_struct scalar A
	A.draws=1		// Will do one draw with no burn as set 
	A.aopt=.234
	A.burn=0
	A.passes=0
	A.totaldraws=0
	A.Accept=0
	A.append="append"
	A.reeval=""
	A.alginfo=amcmc_parser("")
	return(A)
}
void amcmc_lnf(struct amcmc_struct A,pointer (real scalar function) scalar f) A.f=f
void amcmc_args(struct amcmc_struct A,pointer matrix Z)                       A.M=Z
void amcmc_xinit(struct amcmc_struct A,real rowvector xinit) {
																				   A.xinit=xinit
																				   A.xvals=xinit
																				   A.mu=xinit
}
void amcmc_propmean(struct amcmc_struct A,real rowvector mu) A.mu=mu
void amcmc_Vinit(struct amcmc_struct A,real matrix Vinit)                          A.Vinit=Vinit
void amcmc_aopt(struct amcmc_struct A,real scalar aopt)                            A.aopt=aopt
void amcmc_blocks(struct amcmc_struct A,real matrix blocks)                        A.blocks=blocks
void amcmc_model(struct amcmc_struct A,transmorphic M)                             A.M=M
void amcmc_noisy(struct amcmc_struct A,string scalar noisy)                        A.noisy=noisy
void amcmc_alginfo(struct amcmc_struct A,string matrix alginfo)                    A.alginfo=alginfo
void amcmc_damper(struct amcmc_struct A,real scalar damper)                        A.damper=damper
void amcmc_lambda(struct amcmc_struct A,real matrix lam)                           A.lam=lam
void amcmc_draws(struct amcmc_struct A,real scalar draws)                          A.draws=draws
void amcmc_burn(struct amcmc_struct A,real scalar burn)                            A.burn=burn
void amcmc_append(struct amcmc_struct A,string scalar append)                      A.append=append
void amcmc_reeval(struct amcmc_struct A,string scalar reeval)					   A.reeval=reeval
void amcmc_arate(struct amcmc_struct A,real scalar arate)                          A.arate=arate
/* Routines for viewing results and initial settings */
transmorphic amcmc_results_args(struct amcmc_struct A,real scalar i) return((*A.M[i]))
real matrix amcmc_results_draws(struct amcmc_struct A)        return(A.xvals)
real matrix amcmc_results_vals(struct amcmc_struct A)         return(A.vals)
real matrix amcmc_results_arate(struct amcmc_struct A)        return(A.arate)
real matrix amcmc_results_passes(struct amcmc_struct A)       return(A.passes)
real matrix amcmc_results_totaldraws(struct amcmc_struct A)   return(A.totaldraws)
real matrix amcmc_results_aopt(struct amcmc_struct A)         return(A.aopt)
real matrix amcmc_results_acceptances(struct amcmc_struct A)  return(A.Accept)
real matrix amcmc_results_propmean(struct amcmc_struct A)     return(A.mu)
real matrix amcmc_results_propvar(struct amcmc_struct A)      return(A.Vval)
string matrix amcmc_results_alginfo(struct amcmc_struct A)      return(A.alginfo)
string matrix amcmc_results_noisy(struct amcmc_struct A)        return(A.noisy)
real matrix amcmc_results_blocks(struct amcmc_struct A)       return(A.blocks)
real matrix amcmc_results_damper(struct amcmc_struct A)       return(A.damper)
real matrix amcmc_results_xinit(struct amcmc_struct A)        return(A.xinit)
real matrix amcmc_results_Vinit(struct amcmc_struct A)        return(A.Vinit)
real matrix amcmc_results_lambda(struct amcmc_struct A)       return(A.lam)
real matrix amcmc_results_lastdraw(struct amcmc_struct A)     return(A.xvals[rows(A.xvals),])
void amcmc_results_report(struct amcmc_struct A) 
{
printf("Total passes:          %f\n", A.passes)
printf("Iterations per pass:   %f\n", A.draws)
printf("Total iterations:      %f\n",A.totaldraws)
printf("Mean function value:   %f\n",mean(A.vals))
printf("Function val, last it: %f\n",A.vals[rows(A.vals)])
printf("Average accept. rate:  %f\n",rowsum(mean(A.arate)/cols(A.arate)))
printf("Burn-in iterations:    %f\n",A.burn)
}
real scalar amcmc_feval(pointer (real scalar function) scalar f,
					    real rowvector x,
						transmorphic M,
						string scalar etype,
						string scalar opttype)
{
		if      (opttype=="moptimize") return(amcmc_mopt_eval(f,x,M,etype))
		else if (opttype=="optimize")  return(amcmc_opt_eval(f,x,M,etype))
		else if (eltype(M)=="pointer") return(amcmc_eval_args("standalone",f,x,M))	
		else                           return((*f)(x)) 
}
real matrix amcmc_Vsetup(real matrix V,real scalar damper,string scalar fast)
{
	if (fast!="fast" & damper==.) return(cholesky(V)')
	else if (fast=="fast" & rows(V)>1) return(sqrt(diagonal(V))')
	else if (fast=="fast" & rows(V)==1) return(sqrt(V))	
	else return(V)
}
real matrix amcmc_Vdrawprep(real matrix V,real matrix Vold,real scalar damper,string scalar fast)
{
	real matrix Vret
	if (damper==.) Vret=V
	else if (fast=="fast") {
		Vret=sqrt(V)
		if (hasmissing(Vret)) {
			Vret=sqrt(Vold)
			V=Vold
		                    }
						   }
	else {
		Vret=cholesky(V)'
		if (hasmissing(Vret)) {
			Vret=cholesky(Vold)'
			V=Vold
			                }
		 }
	return(Vret)
}
void amcmc_drawinf_update(real matrix mu, real matrix x,real matrix V, real matrix Vold, real matrix lam,real matrix alpha,
                          real scalar aopt, real scalar i,real scalar damper, string scalar fast)
{
	if (damper==.) return
	lam=lam:*exp(1/(i+1)^damper*(alpha:-aopt))
	mu =mu+1/(i+1)^damper*(x-mu)
	Vold=V
	if (fast=="fast") V=V+1/(i+1)^damper*((x-mu):^2-V)	
	else {
		V=V+1/(i+1)^damper*((x-mu)'(x-mu)-V)
		_makesymmetric(V)
         }
}
void amcmc_makenoise(real scalar its,string scalar noisy,real scalar val)
{
	if (round(its/50)==its/50 & noisy=="noisy") {
		printf(" %f: f(x) = %g\n",its,val)
		displayflush()
												}
	else if (noisy=="noisy") {
		printf(".")
		displayflush()
						     }
}
/* Utility functions */
void amcmc_makeargs_1(pointer matrix Z,transmorphic X1) {
	X1=*Z[1]
									   }
void amcmc_makeargs_2(pointer matrix Z,transmorphic X1,
									   transmorphic X2) {
	X1=*Z[1]
	X2=*Z[2]
										   }	
void amcmc_makeargs_3(pointer matrix Z,transmorphic X1,
                                       transmorphic X2,
									   transmorphic X3) {
	X1=*Z[1]
	X2=*Z[2]
	X3=*Z[3]
												 }
void amcmc_makeargs_4(pointer matrix Z,transmorphic X1,
                                       transmorphic X2,
									   transmorphic X3,
									   transmorphic X4) {
	X1=*Z[1]
	X2=*Z[2]
	X3=*Z[3]
	X4=*Z[4]
												 }
void amcmc_makeargs_5(pointer matrix Z,transmorphic X1,
                                       transmorphic X2,
									   transmorphic X3,
									   transmorphic X4,
									   transmorphic X5) {
	X1=*Z[1]
	X2=*Z[2]
	X3=*Z[3]
	X4=*Z[4]
	X5=*Z[5]
												 }
void amcmc_makeargs_6(pointer matrix Z,transmorphic X1,
                                       transmorphic X2,
									   transmorphic X3,
									   transmorphic X4,
									   transmorphic X5,
									   transmorphic X6) {
	X1=*Z[1]
	X2=*Z[2]
	X3=*Z[3]
	X4=*Z[4]
	X5=*Z[5]
	X6=*Z[6]
												 }
void amcmc_makeargs_7(pointer matrix Z,transmorphic X1,
								       transmorphic X2,
									   transmorphic X3,
									   transmorphic X4,
									   transmorphic X5,
									   transmorphic X6,
									   transmorphic X7) {
	X1=*Z[1]
	X2=*Z[2]
	X3=*Z[3]
	X4=*Z[4]
	X5=*Z[5]
	X6=*Z[6]
	X7=*Z[7]
										        }
void amcmc_makeargs_8(pointer matrix Z,transmorphic X1,
                                       transmorphic X2,
									   transmorphic X3,
									   transmorphic X4,
									   transmorphic X5,
									   transmorphic X6,
									   transmorphic X7,
									   transmorphic X8) {
	X1=*Z[1]
	X2=*Z[2]
	X3=*Z[3]
	X4=*Z[4]
	X5=*Z[5]
	X6=*Z[6]
	X7=*Z[7]
	X8=*Z[8]
												}												
void amcmc_makeargs_9(pointer matrix Z,transmorphic X1,
                                       transmorphic X2,
									   transmorphic X3,
									   transmorphic X4,
									   transmorphic X5,
									   transmorphic X6,
									   transmorphic X7,
									   transmorphic X8,
									   transmorphic X9) {
	X1=*Z[1]
	X2=*Z[2]
	X3=*Z[3]
	X4=*Z[4]
	X5=*Z[5]
	X6=*Z[6]
	X7=*Z[7]
	X8=*Z[8]
	X9=*Z[9]
											}						
real scalar amcmc_eval_args(string scalar type,pointer (real scalar function) scalar f,
									 real rowvector b, pointer matrix Z) {
	transmorphic _X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_X7_,_X8_,_X9_,s,H
	real scalar val
	if (length(Z)==1) { 
		amcmc_makeargs_1(Z,_X1_=.)
		if (type=="standalone") return((*f)(b,_X1_))
		(*f)(0,b,_X1_,val=.,s=.,H=.)
		if (type=="d0" | type=="d1" | type=="d2") return(val)
		else {
			if (hasmissing(val)==0) return(sum(val))
			else return(.)
													  }
					}
	else if (length(Z)==2) {
		amcmc_makeargs_2(Z,_X1_=.,_X2_=.)
		if (type=="standalone") return((*f)(b,_X1_,_X2_))
		(*f)(0,b,_X1_,_X2_,val=.,s=.,H=.)
	    if (type=="d0" | type=="d1" | type=="d2") return(val)
		else {
			if (hasmissing(val)==0) return(sum(val))
			else return(.)
													  }
					}
	else if (length(Z)==3) {
		amcmc_makeargs_3(Z,_X1_=.,_X2_=.,_X3_=.)
		if (type=="standalone") return((*f)(b,_X1_,_X2_,_X3_))	
		(*f)(0,b,_X1_,_X2_,_X3_,val=.,s=.,H=.)
		if (type=="d0" | type=="d1" | type=="d2") return(val)
		else {
			if (hasmissing(val)==0) return(sum(val))
			else return(.)
													  }
			       }
	else if (length(Z)==4) {
		amcmc_makeargs_4(Z,_X1_=.,_X2_=.,_X3_=.,_X4_=.)
		if (type=="standalone") return((*f)(b,_X1_,_X2_,_X3_,_X4_))
		(*f)(0,b,_X1_,_X2_,_X3_,_X4_,val=.,s=.,H=.)
        if (type=="d0" | type=="d1" | type=="d2") return(val)
		else {
			if (hasmissing(val)==0) return(sum(val))
			else return(.)
													  }
	                }
	else if (length(Z)==5) {
		amcmc_makeargs_5(Z,_X1_=.,_X2_=.,_X3_=.,_X4_=.,_X5_=.)
		if (type=="standalone") return((*f)(b,_X1_,_X2_,_X3_,_X4_,_X5_))
		(*f)(0,b,_X1_,_X2_,_X3_,_X4_,_X5_,val=.,s=.,H=.)
	    if (type=="d0" | type=="d1" | type=="d2") return(val)
		else {
			if (hasmissing(val)==0) return(sum(val))
			else return(.)
													  }
	                }
	else if (length(Z)==6) {
		amcmc_makeargs_6(Z,_X1_=.,_X2_=.,_X3_=.,_X4_=.,_X5_=.,_X6_=.)
		if (type=="standalone") return((*f)(b,_X1_,_X2_,_X3_,_X4_,_X5_,_X6_))
		(*f)(0,b,_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,val=.,s=.,H=.)
		if (type=="d0" | type=="d1" | type=="d2") return(val)
		else {
			if (hasmissing(val)==0) return(sum(val))
			else return(.)
													  }
	                }
	else if (length(Z)==7) {
		amcmc_makeargs_7(Z,_X1_=.,_X2_=.,_X3_=.,_X4_=.,_X5_=.,_X6_=.,_X7_=.)
		if (type=="standalone") return((*f)(b,_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_X7_))	
		(*f)(0,b,_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_X7_,val=.,s=.,H=.)
		if (type=="d0" | type=="d1" | type=="d2") return(val)
		else {
			if (hasmissing(val)==0) return(sum(val))
			else return(.)
													  }
	                }
	else if (length(Z)==8) {
		amcmc_makeargs_8(Z,_X1_=.,_X2_=.,_X3_=.,_X4_=.,_X5_=.,_X6_=.,_X7_=.,_X8_=.)
		if (type=="standalone") return((*f)(b,_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_X7_,_X8_))	
     	(*f)(0,b,_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_X7_,_X8_,val=.,s=.,H=.)
		if (type=="d0" | type=="d1" | type=="d2") return(val)
		else {
			if (hasmissing(val)==0) return(sum(val))
			else return(.)
													  }
	                } 
	else {
		amcmc_makeargs_9(Z,_X1_=.,_X2_=.,_X3_=.,_X4_=.,_X5_=.,_X6_=.,_X7_=.,_X8_=.,_X9_=.)
		if (type=="standalone") return((*f)(b,_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_X7_,_X8_,_X9_))	
		(*f)(0,b,_X1_,_X2_,_X3_,_X4_,_X5_,_X6_,_X7_,_X8_,_X9_,val=.,s=.,H=.)
        if (type=="d0" | type=="d1" | type=="d2") return(val)
		else {
			if (hasmissing(val)==0) return(sum(val))
			else return(.)
													  }
		}
}
real scalar amcmc_mopt_eval(pointer (real scalar function) scalar f,
							real rowvector beta,
							struct mopt__struct M,
							string scalar evaltype)
{
	real scalar val
	transmorphic s,H
	if (evaltype=="gf0" | evaltype=="gf1" | evaltype=="gf2") {
		(*f)(M,0,beta,val=.,s=.,H=.)
		if (hasmissing(val)==0) return(sum(val))
		else return(.)
															      }
	else if (evaltype=="q") {
	    (*f)(M,0,beta,val=.,s=.,H=.)
		return(val)
						    }
	else {
		(*f)(M,0,beta,val=.,s=.,H=.)
		return(val)
		 }
}
real scalar amcmc_opt_eval(pointer (real scalar function) scalar f,
						   real rowvector beta,
						   struct opt__struct M,
						   string scalar evaltype)
 {
    real scalar nargs,val
	transmorphic Z

	nargs=M.nargs					   
	Z=(M.arglist)[1,1::nargs]
	val=amcmc_eval_args(evaltype,f,beta,Z)
	return(val)
}	
string matrix amcmc_parser(string matrix ainfo)
{
	string scalar type,etype,fast,inine,opttype
	
	if (rowsum(ainfo:=="block")>0)      type="block"
	else if (rowsum(ainfo:=="mwg")>0) type="mwg"
	else type="global"

	if (rowsum(ainfo:=="d0") | rowsum(ainfo:=="d1") | rowsum(ainfo:=="d2")) etype="d0"
	else if (rowsum(ainfo:=="gf0") | rowsum(ainfo:=="gf1") | rowsum(ainfo:=="gf2")) etype="gf0"
	else if (rowsum(ainfo:=="e0") | rowsum(ainfo:=="e1") | rowsum(ainfo:=="e2")) etype="e0"
	else if (rowsum(ainfo:=="v0") | rowsum(ainfo:=="v1") | rowsum(ainfo:=="v2")) etype="v0"	
	else if (rowsum(ainfo:=="q0") | rowsum(ainfo:=="q1") | rowsum(ainfo:=="q2")) etype="q0"
	else etype="standalone"
	
	if (rowsum(ainfo:=="fast")) fast="fast"
	else fast="junk"
	
	if (rowsum(ainfo:=="inline")) inine="inline"
	else inine="junk"
	
	if (rowsum(ainfo:=="optimize")) opttype="optimize"
	else if (rowsum(ainfo:=="moptimize")) opttype="moptimize"
	else opttype="standalone"
	return((type,opttype,etype,fast,inine))
}
mata mlib create lamcmc, dir(PERSONAL) replace
mata mlib add lamcmc *()
mata mlib index
end

