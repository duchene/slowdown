library(phangorn)
library(phytools)
library(geiger)
library(caper)

# The following function is to standardize the values of gamma for a phylogeny by creating random pure birth trees of the same size.

standgam <- function(phy, N = 1000){

	 gams <- vector()
	 ntip <- length(phy$tip.label)
	 for(i in 1:N){
	       gams[i] <- gammaStat(pbtree(n=ntip))
	       
	 }
	 res <- (gammaStat(phy) - mean(gams)) / sd(gams)
	 if(res < qnorm(0.05)){ sig <- 1 } else { sig <- 0 }
	 return(c(res, sig))
}

# The following function receives a single large tree, and extracts a given number of monophyletic clades. These cades have sizes within a given range. Only clades with a given range in latitudes is extracted. The function extracts the species richness, standardized value of gamma, clade age, net diversification rate and mean latitude for each clade.

getLatDat <- function(phy, lats, minsize = 20, maxsize = 200, maxlatrange = 10, N = 100){
	  latdat <- matrix(0, ncol = 5, nrow = 1)
	  colnames(latdat) <- c("lat", "sprich", "zgam", "dr", "age")
	  usedtips <- vector()
	  for(i in 1:N){
	  	newphy <- extract.clade(phy, sample(phy$edge[,1], 1))
		sprich <- length(newphy$tip.label)
                while(sprich < minsize | sprich > maxsize | any(newphy$tip.label %in% usedtips)){
			     newphy <- extract.clade(phy, sample(phy$edge[,1], 1))
                	     sprich <- length(newphy$tip.label)
			     latran <- range(lats[newphy$tip.label])
			     latran <- latran[2] - latran[1]
		}
		#print(sprich)
		usedtips <- append(usedtips, newphy$tip.label[length(newphy$tip.label)])
		phy <- drop.tip(phy, newphy$tip.label[1:(length(newphy$tip.label) - 1)])
		lamfit <- fitContinuous(newphy, lats[newphy$tip.label], model = "lambda")
		rescphy <- rescale(newphy, "lambda", lamfit$opt$lambda )
		lat <- abs(pic(lats[newphy$tip.label], newphy)[1])
		gamres <- standgam(newphy)
		zgam <- gamres[1]
		age <- max(branching.times(newphy))
		dr <- bd.ms(newphy)
		latdat <- rbind(latdat, c(lat, sprich, zgam, dr, age))
		print(paste("done", i,"!"))
	  }
	  latdat <- data.frame(names = usedtips, latdat[2:nrow(latdat),])
	  dat <- list(latdat = latdat, phy = phy)
	  return(dat)

}

# This function uses both functions above to perform the processing of a file with the posterior sample of large trees. For each large tree in the posterior, this function runs PGLS for four primary models. 

getPostreg <- function(postfile, lats, minsize = 20, maxsize = 200, maxlatrange = 10, Ntax = 100, reps = 100){
	   Ntrees <- as.numeric(system(paste("grep -c '.'", postfile), intern = T))
	   mum_mat <- matrix(NA, nrow = Ntrees, ncol = 6)
	   colnames(mum_mat) <- c("mum_sprichlat", "mum_zgamlat", "mum_agelat", "mum_sprichzgam", "mum_sprichage", "mum_zgamage")
	   for(i in 1:Ntrees){
	   	 postre <- read.tree(text = system(paste0("sed '", i,"q;d' ", postfile), intern = T))
	   	 postre <- drop.tip(postre, which(!postre$tip.label %in% names(lats)))
		 m_mat <- matrix(NA, nrow = reps, ncol = 6)
		 colnames(m_mat) <- c("m_sprichlat", "m_zgamlat", "m_agelat", "m_sprichzgam", "m_sprichage", "m_zgamage")
		 for(j in 1:reps){

	   	       regdat <- getLatDat(postre, lats, minsize, maxsize, maxlatrange, Ntax)
		       compdat <- comparative.data(regdat[[2]], regdat[[1]], names.col = "names")
		       res <- try(summary(pgls(sprich ~ lat, compdat, lambda = "ML"))$coef[2, 1], silent = T)
		       if(class(res) == "try-error") print("lm"); res <- summary(lm(sprich ~ lat, data = regdat[[1]]))$coef[2, 1]
		       m_mat[j, 1] <- res
		       res <- try(summary(pgls(zgam ~ lat, compdat, lambda = "ML"))$coef[2, 1], silent = T)
                       if(class(res) == "try-error") print("lm"); res <- summary(lm(zgam ~ lat, data = regdat[[1]]))$coef[2, 1]
                       m_mat[j, 2] <- res
		       res <- try(summary(pgls(age ~ lat, compdat, lambda = "ML"))$coef[2, 1], silent = T)
                       if(class(res) == "try-error") print("lm"); res <- summary(lm(age ~ lat, data = regdat[[1]]))$coef[2, 1]
                       m_mat[j, 3] <- res
		       res <- try(summary(pgls(sprich ~ zgam, compdat, lambda = "ML"))$coef[2, 1], silent = T)
                       if(class(res) == "try-error") print("lm"); res <- summary(lm(sprich ~ zgam, data = regdat[[1]]))$coef[2, 1]
                       m_mat[j, 4] <- res
		       res <- try(summary(pgls(sprich ~ age, compdat, lambda = "ML"))$coef[2, 1], silent = T)
                       if(class(res) == "try-error") print("lm"); res <- summary(lm(sprich ~ age, data = regdat[[1]]))$coef[2, 1]
                       m_mat[j, 5] <- res
		       res <- try(summary(pgls(zgam ~ age, compdat, lambda = "ML"))$coef[2, 1], silent = T)
                       if(class(res) == "try-error") print("lm"); res <- summary(lm(zgam ~ age, data = regdat[[1]]))$coef[2, 1]
                       m_mat[j, 6] <- res

	   	 }
		 mum_mat[i,] <- colMeans(m_mat)
		 print(paste("Posterior tree", i, "processed"))
		 
	   }

	   return(mum_mat)

}