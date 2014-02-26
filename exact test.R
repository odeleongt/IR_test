M <- 114
N <- 2942
alpha <- 0.05
setwd("D:	//uvg//0_Tesis//Rfun")

ns <- 
{
	ns <- round(rbeta(M,0.09,0.12)*49+1)
	
	sumns <- sum(ns)

	while (sumns!=N)
	{
		dif <- sumns-N
		signflag <- sign(dif)
		dif <- abs(dif)
		viable <- which(ns>1 & ns<50)
		if(length(viable)<1) viable <- 1:M
		mod <- round(runif(dif,1,length(viable)))
		ns[viable[mod]] <- ns[viable[mod]] - signflag
		ns[ns>50] <- 50
		ns[ns<1] <- 1
		sumns <- sum(ns)
	}
	ns
}

mns <- hist(ns, plot=FALSE, 0:50)$counts
mni <- rep(1:50, mns)

po <- 1/10000
pa <- 1/1000

k <- sum(mns>0)

bs <- function(N, M, p)
{
	ixs <- rbinom(N,1,p)
	re <- rep(1:M, ns)
	xs <- rep(0,M)
	xs[re[ixs==1]] <- 1
	
	xs
}
#--- Fin de bs ------------------------


X2a1 <- qchisq(alpha, 1) # Valor crítico de chi-cuadrado, nivel alpha, 1 grado de libertad


dlogLdp <- function(p,x)
{
	p1 <- 1-p
	sumi <- 0
	ks <- which(mns>0)
	mis <- mns[mns>0]
	for(i in 1:k)
	{
		ni <- ks[i]
		xis <- x[mni==ni]
		sumj <- 0
		j <- mis[i]
		while(j>0)
		{
			xij <- xis[j]
			sumj <- sumj + (((xij*ni*p1^(ni-1))/(1-p1^ni))-((ni*(1-xij))/(p1)))
			j <- j-1
		}
		sumi <- sumi + sumj
	}
	sumi
}


# Fisher's information
I <- function(p)
{
	i <- 1:k
	sum(mns[i]*(((i^2)*((1-p)^(i-2)))/(1-(1-p)^i)))
}


scoretest <- function(po,x)
{
	(dlogLdp(po,x)/sqrt(I(po)))^2
}



dos <- function(probis, times)
{
stime <- array(dim=c(49+7,3))
ins <- 1:length(probis)

for(souterindex in ins)
{
a <- Sys.time()
po <- probis[souterindex]


#### Simulación 1 ---------------------
scorespo <- NULL

for(index in 1:times)
{
	x <- bs(N,M,po)
	if(sum(x)==0) x[round(runif(1,1,M))] <- 1
	scorespo[index] <- scoretest(po,x)
}
while(sum(is.na(scorespo))>0)
{
	for(i in which(is.na(scorespo)==TRUE))
	{
		x <- bs(N,M,po)
		if(sum(x)==0) x[round(runif(1,1,M))] <- 1
		scorespo[i] <- scoretest(po,x)
	}
}
write.csv(scorespo,paste("p=",po,"//Scores//scores, po=",po,".csv", sep=""))
scoresqs <- quantile(scorespo,probs= c(0.025,0.975), na.rm=TRUE)
write.csv(scoresqs,paste("p=",po,"//Scores//quantiles, po=",po,".txt", sep=""))
b <- Sys.time()
c <- b-a
stime[((souterindex-1)*8)+1,] <- c("-----","-----","---------------")
stime[((souterindex-1)*8)+souterindex+1,] <- c(po,"--",c)

for(sinnerindex in ins[-souterindex])
{
pa <- probis[sinnerindex]
a <- Sys.time()

#### Simulación 2 ---------------------
scorespa <- NULL
rejects <- NULL

for(i in 1:times)
{
	x <- bs(N,M,pa)
	if(sum(x)==0) x[round(runif(1,1,M))] <- 1
	scorespa[i] <- scoretest(po,x)
}

while(sum(is.na(scorespa))>0)
{
	for(i in which(is.na(scorespa)==TRUE))
	{
		x <- bs(N,M,pa)
		if(sum(x)==0) x[round(runif(1,1,M))] <- 1
		scorespa[i] <- scoretest(po,x)
	}
}

rejects <- scorespa>X2a1
write.csv(scorespa,paste("p=",po,"//Scores//scores, pa=",pa,".csv", sep=""))

spower <- (sum(rejects)*100)/times

write.csv(spower,paste("p=",po,"//Scores//po=",po,", pa=",pa,", power.txt", sep=""))

b <- Sys.time()
c <- b-a
stime[((souterindex-1)*8)+sinnerindex+1,] <- c(po,pa,c)
}
}
stime
}

ps <- c(1/10000, 1/1000)
times <- 100
tiempos <- dos(probis=ps, times=times)
write.csv(tiempos,"scoresTimes.csv")
tiempos
