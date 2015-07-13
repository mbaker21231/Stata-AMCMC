use http://fmwww.bc.edu/repec/bocode/t/traindata.dta
gen mprice = -price
mixlogitwtp y contract local wknown, group(gid) id(pid) price(mprice) ///
            rand(tod seasonal) nrep(500)
bayesmixedlogitwtp y contract local wknown, group(gid) id(pid) price(price) ///
		    rand(tod seasonal) noisy draws(4000) burn(2000) thin(4) jumble
bayesmixedlogitwtp y, group(gid) id(pid) price(price) ///
		    rand(tod seasonal wknown) noisy draws(500) burn(250) 
bayesmixedlogitwtp y, group(gid) id(pid) price(price) ///
		    rand(seasonal tod) samplerr("mwg") noisy draws(500) burn(250) 
mixlogitwtp y, group(gid) id(pid) price(mprice) ///
            rand(seasonal tod) nrep(500)

			/* Seem to compare well */
			