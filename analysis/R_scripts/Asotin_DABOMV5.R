
### FIX CODE IN LINE 140 and init_fnc for AFC####
rm(list=ls())
ls()


library(coda)
library(R2jags)
# library(runjags)
# library(rjags)
library(mcmcplots)


setwd("C:\\data\\Asotin_Code\\AdultCJS\\Papers\\AbundanceDiversity\\") # change this to your directory

#this is and individual CH for all sites in 2021 for wild females released at ASOTIC
#low water year detection is near 100%, added dummy data for CCA and George since we had no detections to test code
#this version uses 3 detections in CCA
#this version uses AFC mainstem, NF, SF

dat<- read.csv("Asotin2021FemalesDABOM.csv",header=TRUE)
dat<-dat[,-1]
datfinal<-dat
colnames(datfinal)<-NULL


#data
#matrix with Geo, ACBBO, ACBA0, CCABO, CCAMO, CCAAO, AFCBO, AFFNO, AFFSO
cap_hist <- as.matrix(datfinal)
n_fish<-dim(cap_hist)[1]
n_branch_AFC =3
n_branch_ACB =3
n_branch_ASOTIC=3
zero_vec<-rep(0, max(n_branch_ACB, n_branch_ASOTIC,n_branch_AFC))
ACB_dirch_vec = rep(0.01, n_branch_ACB) #changed from 1 to 0.01
ASOTIC_dirch_vec = rep(0.01, n_branch_ASOTIC)
AFC_dirch_vec = rep(0.01, n_branch_ACB) #changed from 1 to 0.01

#beta priors for capture probs changed from 1
a<- 0.01
b<-0.01

#error message Unknown Variable ACB_dirch_vec
#likely error message Unknown Variable ASOTic_dirch_vec

jags.data <- list(cap_hist=cap_hist,
                  n_fish=n_fish,
                  n_branch_AFC =n_branch_AFC,
                  n_branch_ACB =n_branch_ACB,
                  n_branch_ASOTIC=n_branch_ASOTIC,
                  AFC_dirch_vec = AFC_dirch_vec,
                  ACB_dirch_vec = ACB_dirch_vec,
                  ASOTIC_dirch_vec = ASOTIC_dirch_vec,
                  zero_vec=zero_vec,
                  a=a,b=b)#priors for detection 

# set initial values to avoid "Invalid parent values" errors
init_fnc = function() {
  list(a_ASOTIC = ifelse(cap_hist[,1] == 1,
                                    1,
                                   ifelse(rowSums(cap_hist[,c(2:9)]) > 1,
                                          2, 3)),
                 a_ACB = ifelse(rowSums(cap_hist[,c(4,5,6)]) > 1,
                                1,
                                ifelse(rowSums(cap_hist[,c(7:9)]) > 1,
                                       2, 3)),
                a_AFC = ifelse(rowSums(cap_hist[,c(7,8)]) > 1,
                                1,
                                ifelse(rowSums(cap_hist[,c(8,9)]) > 1,
                                      2, 3)))

  }
#init_fnc()


sink("AsotinDABOM.jags")
cat("
model {
  # Code for Kevin See
  #psi[f,l] = occurence probability that fish f actually passed location l
  #p = probability of detection/capture for each site
  #psi_Asotin[1]=trans to GEORC
  #psi_Asotin[2]=trans to ACB
  #psi_Asotin[3]=last obs above ASOTIC and below GEORGC and below ACB
  #psi_ACB[1]=trans to CCA
  #psi_ACB[2]=trans to ACB
  #psi_ACB[3]=last obs above ACB and below CCA and below AFC
  #psi_AFC[1]=trans to AFC South
  #psi_AFC[2]=trans to AFC North
  #psi_AFC[3]=last obs above AFC mainstem and below AFC north and south
 
  

  # Priors for detection probabilities
  GEORGC_p <- 1 # assume perfect detection
  ACBB0_p ~ dbeta(a,b)T(0.0001 ,0.9999);
  ACBA0_p ~ dbeta(a,b)T(0.0001 ,0.9999);
  CCAB0_p ~ dbeta(a,b)T(0.0001 ,0.9999);
  CCAM0_p ~ dbeta(a,b)T(0.0001 ,0.9999);
  CCAA0_p ~ dbeta(a,b)T(0.0001 ,0.9999);
  AFCB0_p ~ dbeta(a,b)T(0.0001 ,0.9999);
  AFCS0_p <-1
  AFCN0_p<- 1


  # Priors for transition probabilities
  psi_ASOTIC[1:n_branch_ASOTIC] ~ ddirch(ASOTIC_dirch_vec);
  omega_ASOTIC[1, 1:n_branch_ASOTIC] <- zero_vec[1:(n_branch_ASOTIC)];
  omega_ASOTIC[1, (n_branch_ASOTIC + 1)] <- 1;
  omega_ASOTIC[2, 1:n_branch_ASOTIC] <- psi_ASOTIC;
  omega_ASOTIC[2, (n_branch_ASOTIC + 1)] <- 0;

  psi_ACB[1:n_branch_ACB] ~ ddirch(ACB_dirch_vec);
  omega_ACB[1, 1:n_branch_ACB] <- zero_vec[1:(n_branch_ACB)]; # fish that didn't go past ACB
  omega_ACB[1, (n_branch_ACB + 1)] <- 1;
  omega_ACB[2, 1:n_branch_ACB] <- psi_ACB;
  omega_ACB[2, (n_branch_ACB + 1)] <- 0;
  
  #check these
  psi_AFC[1:n_branch_AFC] ~ ddirch(AFC_dirch_vec);
  omega_AFC[1, 1:n_branch_AFC] <- zero_vec[1:(n_branch_AFC)]; 
  omega_AFC[1, (n_branch_AFC + 1)] <- 1;
  omega_AFC[2, 1:n_branch_AFC] <- psi_AFC;
  omega_AFC[2, (n_branch_AFC + 1)] <- 0;
  

  # Where is each fish?
  for(i in 1:n_fish) {
    a_ASOTIC[i] ~ dcat( psi_ASOTIC );
    for (j in 1:n_branch_ASOTIC)	{
      eta_ASOTIC[i,j] <- equals(a_ASOTIC[i],j) # equals(x,y) is a test for equality, returns [1,0]
    }

    # the row number acts as switch between rows 1&2 using stochastic node
    a_ACB[i] ~ dcat( omega_ACB[(eta_ASOTIC[i,2] + 1), 1:(n_branch_ACB+1)] )
    for (j in 1:n_branch_ACB)	{
      eta_ACB[i,j] <- equals(a_ACB[i],j) 
    }
 
    #Check this
    a_AFC[i] ~ dcat( omega_AFC[(eta_ACB[i,2] + 1), 1:(n_branch_AFC+1)] )
    for (j in 1:n_branch_AFC)	{
      eta_AFC[i,j] <- equals(a_AFC[i],j) 
    }
  } # end the n_fish loop


  # Were tags observed?
  # Check cap_hist[,6:9]
  for (i in 1:n_fish) {
    cap_hist[i,1] ~ dbern( GEORGC_p * eta_ASOTIC[i,1] );
    cap_hist[i,2] ~ dbern( ACBB0_p * eta_ASOTIC[i,2] );
    cap_hist[i,3] ~ dbern( ACBA0_p * eta_ASOTIC[i,2] );
    cap_hist[i,4] ~ dbern( CCAB0_p * eta_ACB[i,1] );
    cap_hist[i,5] ~ dbern( CCAM0_p * eta_ACB[i,1] );   
    cap_hist[i,6] ~ dbern( CCAA0_p * eta_ACB[i,1] );
    cap_hist[i,7] ~ dbern( AFCB0_p * eta_ACB[i,2] );
    cap_hist[i,8] ~ dbern( AFCS0_p * eta_AFC[i,1] );
    cap_hist[i,9] ~ dbern( AFCN0_p * eta_AFC[i,2] );    
  }  # end the n_fish loop

 #Calculate last location here

}
",fill = TRUE)
sink()

#removed parallel so this should run in R2jags
start = Sys.time()
out<- jags(data=jags.data, model.file="AsotinDABOM.jags",inits=init_fnc,
                    n.chains=4, n.thin=10,n.burnin=10000,n.iter=51000,
                    parameters.to.save=c( "GEORGC_p","ACBB0_p","ACBA0_p","CCAA0_p","CCAM0_p","CCAB0_p","AFCA0_p","AFCB0_p",
                                         "psi_ASOTIC","psi_ACB") )
(elapsed = Sys.time() - start)

print(out,dig=3)




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



summary(out)

# next step is to multiply some movement parameters together appropriately
post_mat = as.matrix(out,
                     iters = T,
                     chains = T)

# pull out relevant columns (transition probabilities for wild fish)
trans_mat <- as.data.frame(post_mat[,c(1,2,
                                       13:15,
                                       10:12)])

# multiply the transistion probabilities above ACB by the probability of a fish getting to ACB in the first place (psi_ASOTIC[2])
trans_mat$`psi_ACB[1]` = trans_mat$`psi_ACB[1]` * trans_mat$`psi_ASOTIC[2]`
trans_mat$`psi_ACB[2]` = trans_mat$`psi_ACB[2]` * trans_mat$`psi_ASOTIC[2]`
trans_mat$`psi_ACB[3]` = trans_mat$`psi_ACB[3]` * trans_mat$`psi_ASOTIC[2]`

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
summary(mcmc_lst)
head(mcmc_lst)