require(rjags)
require(coda)

setwd("C:/Users/ndhen/Dropbox (UW)/School/Misc/JAGS tutorial/Joseph - NGS/jags_tutorial_joseph/")

# simulate data
stool    <- c(rep(1, 40), rep(0, 122))
serology <- c(rep(1, 38), rep(0, 2), rep(1, 87), rep(0, 35))
N        <- 162
set.seed(108)

# create the model for stool exam
model_string_stool = "
model {
  for(i in 1:N){
    y[i] ~ dbern(S*p+(1-C)*(1-p))
  }
  S ~ dbeta(4.44, 13.31)
  C ~ dbeta(71.25, 3.75)
  p ~ dunif(0,1)
}
"

# create the model for serology
model_string_sero = "
model {
  for(i in 1:N){
    y[i] ~ dbern(S*p+(1-C)*(1-p))
  }
  S ~ dbeta(21.96, 5.49)
  C ~ dbeta(4.1, 1.76)
  p ~ dunif(0, 1)
}
"

# write model to temporary file
writeLines(model_string_stool, con = "stool_model.txt")
writeLines(model_string_sero, con  = "sero_model.txt")

# initialize the model 
jags_stool <- jags.model("stool_model.txt",
                         data     = list('y' = stool,
                                         'N' = N),
                         n.chains = 4,
                         n.adapt  = 1000)
update(jags_stool, 10000) #burn in

stool_samples <- coda.samples(jags_stool, 
                              c("p", "S", "C"), 
                              n.iter = 10000)
summary(stool_samples)

##############

# initialize the model 
jags_sero <- jags.model("sero_model.txt",
                         data     = list('y' = serology,
                                         'N' = N),
                         n.chains = 4,
                         n.adapt  = 1000)
update(jags_sero, 10000) #burn in

sero_samples <- coda.samples(jags_sero, 
                             c("p", "S", "C"), 
                             n.iter = 10000)
summary(sero_samples)

# combine plots
stool_samples_df <- as.data.frame(as.matrix(stool_samples))
sero_samples_df  <- as.data.frame(as.matrix(sero_samples))
plot(density(stool_samples_df$S),
     xlim = c(0, 1),
     ylim = c(0, 14),
     xlab = "Sensitivities and Specificities",
     ylab = "Posterior density",
     lty  = 1,
     main = "")
lines(density(sero_samples_df$S),
      lty = 2)
lines(density(stool_samples_df$C),
      lty = 3)
lines(density(sero_samples_df$C),
      lty = 4)
legend(0, 14,
       legend = c("Sensitivity of stool examination",
                  "Specificity of stool examination",
                  "Sensitivity of serologic test",
                  "Specificity of serologic test"),
       lty=1:4,
       cex=0.8)

####################
# Two-test version #
####################

ones  <- rep(1, times = N)
tests <- data.frame("stool"    = stool,
                    "serology" = serology)

model_string_two_test <- "
var pr[N], q[N,4]
model {
  for(i in 1:N){
    q[i, 1] <- p*(S1*S2) + (1-p)*((1-C1)*(1-C2)) #Pr(+,+)
    q[i, 2] <- p*(S1*(1-S2)) + (1-p)*((1-C1)*C2) #Pr(+,-)
    q[i, 3] <- p*((1-S1)*S2) + (1-p)*(C1*(1-C2)) #Pr(-,+)
    q[i, 4] <- p*((1-S1)*(1-S2)) + (1-p)*(C1*C2) #Pr(-,-)

  L[i] <- equals(tests[i,1],1) * equals(tests[i,2],1) * q[i,1]
        + equals(tests[i,1],1) * equals(tests[i,2],0) * q[i,2]
        + equals(tests[i,1],0) * equals(tests[i,2],1) * q[i,3]
        + equals(tests[i,1],0) * equals(tests[i,2],0) * q[i,4]
  
  pr[i] <- L[i] / 1
  ones[i] ~ dbern(pr[i])
  }

  S1 ~ dbeta(4.44, 13.31)
  C1 ~ dbeta(71.25, 3.75)
  S2 ~ dbeta(21.96, 5.49)
  C2 ~ dbeta(4.1, 1.76)
  p ~ dunif(0,1)
}
"

writeLines(model_string_two_test, con = "two_test_model.txt")
jags_two_test <- jags.model("two_test_model.txt",
                        data     = list('tests' = tests,
                                        'ones'  = ones,
                                        'N'     = N),
                        n.chains = 4,
                        n.adapt  = 1000)
update(jags_two_test, 10000) #burn in

two_test_samples <- coda.samples(jags_two_test, 
                             c("p", "S1", "C1", "S2", "C2"), 
                             n.iter = 10000)
summary(two_test_samples)

two_test_samples_df  <- as.data.frame(as.matrix(two_test_samples))
plot(density(stool_samples_df$p),
     xlim = c(0, 1),
     ylim = c(0, 5),
     xlab = "Prevalence",
     ylab = "Posterior density",
     lty  = 1,
     main = "")
lines(density(sero_samples_df$p),
      lty = 2)
lines(density(two_test_samples_df$p),
      lty = 3)
abline(h  = 1,
      lty = 4)
legend(0, 4.5,
       legend = c("Stool examination",
                  "Serologic test",
                  "Both tests combined",
                  "Prior"),
       lty=1:4,
       cex=0.8)