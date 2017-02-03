# Relative to low-richness habitats, are observations of colonizations and extinctions more likely in high-richness habitats, even after linearly controlling for richness?
# Does the probability of observing (p_o()) colonization (c) or extinction (e) events vary (correlate) with richness (R)? (1) c ~ R?  (2) e ~? (3) c/R ~ R? (4) e/R ~ R?
# 




# ============
# = Gompertz =
# ============
gomp_mean <- function(a, c=b+1, b=c-1){a/(1-c)}
gomp_var <- function(a, sigma, c=b+1, b=c-1){sigma^2/(1-c^2)}

gomp_a <- function(mu, c=b+1, b=c-1){mu*(1-c)}
gomp_c <- function(mu, a){1-a/mu}

gomp_delt <- function(N, a, c=b+1, b=c-1, sigma, log=FALSE){
	n <- length(N)
	PE <- rnorm(n=n, mean=0, sd=sigma)
	if(log){
		a + (b+1)*N + PE
	}else{
		N*exp(a + b*log(N) + PE)
	}
}


# ===============
# = Observation =
# ===============
veil <- function(prob, effort=20){
	#' @param prob probaility of observing a species
	
}


# ==========================
# = Dynamic Classification =
# ==========================
is.ext <- function(x){diff(x>0)<0}
is.col <- function(x){diff(x>0)>0}


# ===============================================================================
# = Extinction Rate 1: Abundance is log-normal, extinction rate is logit(abund) =
# ===============================================================================
eR1 <- function(R, reps=5E2){
	mean( # take the mean of the replicates
		colSums( # sum number of extinctions across the 10 spp
			replicate( # replicate abund & ext realization (p(ext|abund) is deterministic)
				reps, rbinom( # extinctions realized as bernoulli 
					R, 1, 1-plogis( # extinction probability a function of abundance, via 1-plogis(abundance)
						rlnorm( # species abundance is log-normally distributed
							R), 0)))))
}

Richness <- 10:200
extinctionRate1 <- sapply(Richness, eR1)
plot(Richness, extinctionRate1)
summary(lm(extinctionRate1~Richness))


# ==================================================
# = Extinction Rate 1.5: Quick Alternative with Ed =
# ==================================================
# extinction probability declines linearly with abundance
alt <- function(x) x/max(x)
eR1.5 <- function(R, reps=5E2){
	mean( # take the mean of the replicates
		colSums( # sum number of extinctions across the 10 spp
			replicate( # replicate abund & ext realization (p(ext|abund) is deterministic)
				reps, rbinom( # extinctions realized as bernoulli 
					R, 1, 1 - alt( # extinction probability a function of abundance, via 1-alt(abundance)
						rlnorm( # species abundance is log-normally distributed
							R))))))
}

Richness <- 10:200
extinctionRate1.5 <- sapply(Richness, eR1.5)
plot(Richness, extinctionRate1.5)
summary(lm(extinctionRate1.5~Richness))


# =======================================================================
# = Extinction Rate 2: Abund~logN, eR = exp(-lambda*t), lambda=f(b,n,k) =
# =======================================================================
# Model on the relationship between abundance and detectability (extinction) taken from:
# McCarthy et al. The influence of abundance on detectability (2013) Oikos 122:717-726, doi:10.1111/j.1600-0706.2012.20781.x
logLambda <- function(b, abund, k){b*log(abund) + log(k)}
eR2 <- function(R, b, abund, k, dt=1, ..., reps=5E2){
	if(missing(b)){b <- 1} #c(0.5, 1, 1.5)}
	if(missing(abund)){abund <- rlnorm(n=R)}
	if(missing(k)){k <- 0.5}#rep(0.5, length(b))}
	
	lambda <- exp(logLambda(b=b, abund=abund, k=k))
	prExt <- exp(-lambda*dt)
	
	return(sum(prExt))
}

params_eR2 <- as.data.table(expand.grid(R=(1:500), b=c(0.15, 1, 1.85), k=0.5))
params_eR2[,nExt:=mean(replicate(1E2,eR2(R=R, b=b,k=k))),by=c(names(params_eR2))]
# params_eR2[,image(z=matrix(nExt, ncol=length(unique(b))), x=unique(R), y=unique(b))] # hard to determine slopes
params_eR2[,j={
	ylim <- range(nExt)
	u_b <- unique(b)
	M <- length(u_b)
	b_cols <- rainbow(n=M) #magma(M)
	
	.SD[b==u_b[1],j={plot(R, nExt, type='l', ylim=ylim, col=b_cols[1])}]
	for(m in 2:M){
			.SD[b==u_b[m],j={lines(R, nExt, col=b_cols[m])}]
	}
}]







