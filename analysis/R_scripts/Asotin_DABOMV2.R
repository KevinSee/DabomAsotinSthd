

rm(list=ls())
ls()


library(coda)
library(R2jags)
# library(runjags)
library(rjags)
library(mcmcplots)


setwd("C:\\data\\Asotin_Code\\AdultCJS\\Papers\\Abundance and Diversity\\") # change this to your directory

#this is and individual CH for all sites in 2021 for wild females released at ASOTIC
#low water year detection is near 100%, added dummy data for CCA and George since we had no detections to test code

dat<- read.csv("Asotin2021FemalesDABOM.csv",header=TRUE)
# dat<- read.csv("data/Asotin2021FemalesDABOM.csv",header=TRUE)
#head(dat)


#combine AFC/north and AFC/South detections
dat[,10]<-dat[,8]+dat[,9]
dat10<-replace(dat[,10],dat[,10]>1,1)
dat<-dat[,-c(8,9,10)]
dat<-cbind(dat,dat10)

#combine CCA.Down and CCA.Mid
dat[,9]<-dat[,4]+dat[,5]
datCCA<- replace(dat[,9],dat[,9]>1,1)
dat1<-dat[,c(1:3)]
dat2<-dat[,c(6:8)]
datfinal<-cbind(dat1,datCCA)
datfinal<-cbind(datfinal,dat2)
colnames(datfinal)<-NULL

#data
#matrix with Geo,ACBBO,ACA),CCABO,CCAAO,AFCBO,AFFAO
cap_hist <- as.matrix(datfinal)
n_fish<-dim(cap_hist)[1]
fish_type= rep(1, n_fish) # all wlld females
n_branch_ACB =3
n_branch_ASOTIC=3
zero_vec<-rep(0, max(n_branch_ACB, n_branch_ASOTIC))
ACB_dirch_vec = matrix(c(rep(1, n_branch_ACB),
                         rep(1, n_branch_ACB)),
                         nrow = 2,
                       byrow = T)
ASOTIC_dirch_vec = matrix(c(rep(1, n_branch_ASOTIC),
                            rep(1, n_branch_ASOTIC)),
                          nrow = 2,
                          byrow = T)

#error message Unknown Variable ACB_dirch_vec
#likely error message Unknown Variable ASOTic_dirch_vec

jags.data <- list(cap_hist=cap_hist,
                  n_fish=n_fish,
                  fish_type=fish_type,
                  n_branch_ACB =n_branch_ACB,
                  n_branch_ASOTIC=n_branch_ASOTIC,
                  ACB_dirch_vec = ACB_dirch_vec,
                  ASOTIC_dirch_vec = ASOTIC_dirch_vec,
                  zero_vec=zero_vec)

# set initial values to avoid "Invalid parent values" errors
init_fnc = function() {
  list(a_ASOTIC = ifelse(cap_hist[,1] == 1,
                                    1,
                                   ifelse(rowSums(cap_hist[,c(2:7)]) > 1,
                                          2, 3)),
                 a_ACB = ifelse(rowSums(cap_hist[,c(4,5)]) > 1,
                                1,
                                ifelse(rowSums(cap_hist[,c(6:7)]) > 1,
                                       2, 3)))
}


sink("AsotinDABOM.jags")
cat("
model {
  # Code for Kevin See
  #psi[f,l] = occurence probability that fish f actually passed location l
  #p = probability of detection/capture for each site
  #omega = branches upstream

  # Priors for detection probabilities
  GEORGC_p <- 1 # assume perfect detection
  ACBB0_p ~ dbeta(1, 1);
  ACBA0_p ~ dbeta(1, 1);
  CCAB0_p ~ dbeta(1, 1);
  CCAA0_p ~ dbeta(1, 1);
  AFCB0_p ~ dbeta(1, 1);
  AFCA0_p ~ dbeta(1, 1);

  # Priors for transition probabilities
  psi_ASOTIC[1, 1:n_branch_ASOTIC] ~ ddirch(ASOTIC_dirch_vec[1,]);
  psi_ASOTIC[2, 1:n_branch_ASOTIC] ~ ddirch(ASOTIC_dirch_vec[2,]);

  omega_ASOTIC[1, 1:n_branch_ASOTIC] <- zero_vec[1:(n_branch_ASOTIC)];
  omega_ASOTIC[1, (n_branch_ASOTIC + 1)] <- 1;

  omega_ASOTIC[2, 1:n_branch_ASOTIC] <- psi_ASOTIC[1,];
  omega_ASOTIC[2, (n_branch_ASOTIC + 1)] <- 0;

  omega_ASOTIC[3, 1:n_branch_ASOTIC] <- psi_ASOTIC[2,];
  omega_ASOTIC[3, (n_branch_ASOTIC + 1)] <- 0;


  psi_ACB[1, 1:n_branch_ACB] ~ ddirch(ACB_dirch_vec[1,]);
  psi_ACB[2, 1:n_branch_ACB] ~ ddirch(ACB_dirch_vec[2,]);

  omega_ACB[1, 1:n_branch_ACB] <- zero_vec[1:(n_branch_ACB)];
  omega_ACB[1, (n_branch_ACB + 1)] <- 1;

  omega_ACB[2, 1:n_branch_ACB] <- psi_ACB[1,];
  omega_ACB[2, (n_branch_ACB + 1)] <- 0;

  omega_ACB[3, 1:n_branch_ACB] <- psi_ACB[2,];
  omega_ACB[3, (n_branch_ACB + 1)] <- 0;

  # Where is each fish?
  for(i in 1:n_fish) {
    a_ASOTIC[i] ~ dcat( omega_ASOTIC[fish_type[i] + 1, 1:(n_branch_ASOTIC+1)] )
    for (j in 1:n_branch_ASOTIC)	{
      eta_ASOTIC[i,j] <- equals(a_ASOTIC[i],j) # equals(x,y) is a test for equality, returns [1,0]
    }

    a_ACB[i] ~ dcat( omega_ACB[(eta_ASOTIC[i,2] * fish_type[i] + 1), 1:(n_branch_ACB+1)] )
    for (j in 1:n_branch_ACB)	{
      eta_ACB[i,j] <- equals(a_ACB[i],j) # equals(x,y) is a test for equality, returns [1,0]
    }
  } # end the n_fish loop


  # Were tags observed?
  for (i in 1:n_fish) {
    cap_hist[i,1] ~ dbern( GEORGC_p * eta_ASOTIC[i,1] );
    cap_hist[i,2] ~ dbern( ACBB0_p * eta_ASOTIC[i,2] );
    cap_hist[i,3] ~ dbern( ACBA0_p * eta_ASOTIC[i,2] );
    cap_hist[i,4] ~ dbern( CCAB0_p * eta_ACB[i,1] );
    cap_hist[i,5] ~ dbern( CCAA0_p * eta_ACB[i,1] );
    cap_hist[i,6] ~ dbern( AFCB0_p * eta_ACB[i,2] );
    cap_hist[i,7] ~ dbern( AFCA0_p * eta_ACB[i,2] );
  }  # end the n_fish loop
}
",fill = TRUE)
sink()

start.time<-Sys.time()
out<- jags.parallel(data=jags.data, model.file="AsotinDABOM.jags",#inits=inits,
                    n.chains=4, n.thin=1,n.burnin=1000,n.iter=5000,
                    parameters.to.save=c( "GEORGC_p","ACBB0_p","ACBA0_p","CCAA0_p","CCAB0_p","AFCA0_p","AFCB0_p",
                                         "psi_ASOTIC","psi_ACB" ) )#"cp"

# run the model
# i had trouble installing R2jags, so used rjags instead
jags <- jags.model("AsotinDABOM.jags",
                  data = jags.data,
                  inits = init_fnc,
                  n.chains = 4,
                  n.adapt = 5000)
out <- coda.samples(jags,
                    c( "GEORGC_p","ACBB0_p","ACBA0_p","CCAA0_p","CCAB0_p","AFCA0_p","AFCB0_p",
                       "psi_ASOTIC","psi_ACB" ),
                    n.iter = 5000,
                    thin = 5)


end.time<-Sys.time()
model.time<-end.time-start.time
print(model.time)

print(out,dig=3)

# next step is to multiply some movement parameters together appropriately
post_mat = as.matrix(out,
                     iters = T,
                     chains = T)

# pull out relevant columns (transition probabilities for wild fish)
trans_mat <- as.data.frame(post_mat[,c(1,2,
                                       16, 18, 20,
                                       10, 12, 14)])

# multiply the transistion probabilities above ACB by the probability of a fish getting to ACB in the first place (psi_ASOTIC[1,2])
trans_mat$`psi_ACB[1,1]` = trans_mat$`psi_ACB[1,1]` * trans_mat$`psi_ASOTIC[1,2]`
trans_mat$`psi_ACB[1,2]` = trans_mat$`psi_ACB[1,2]` * trans_mat$`psi_ASOTIC[1,2]`
trans_mat$`psi_ACB[1,3]` = trans_mat$`psi_ACB[1,3]` * trans_mat$`psi_ASOTIC[1,2]`

# rename the columns to reflect what point they represent movement past
colnames(trans_mat) = c("CHAIN", "ITER",
                        "GEORGC",
                        "ACB",
                        "ASOTIC_bb",
                        "CCA",
                        "AFC",
                        "ACB_bb")

# put it back into MCMC list if that's easier (including detection parameters)
mcmc_lst = as.mcmc.list(lapply(split(cbind(post_mat[,c(3:9)], trans_mat[,-c(1,2)]),
                                     trans_mat$CHAIN),
                               mcmc))
