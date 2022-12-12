## Author : Nicholas Siame Adam  					          ########
## Code name: Data simulation for MSc Biostatistics ########
## Last update date: 11.11.2021 			          		########
############################################################

library(rms) #regression modeling strategies(works in missing data, non-linear transformations)
library(lmerTest) #new modeling in leaner mixed models (fits with lme4) with p-values
library(pbkrtest) #new Methods for Tests in Linear Mixed Models
library(Metrics) #implements metrics for regression or information retrieval
install.packages("hmer") 
library(hmer)
(options(timeout = 600)) # increasing package installation timeout 

# Generating objects to store the results

#for clusters, response and explanatory variables
clusterId<-c()
Y<-c()
X1<-c()
X2<-c()
U<-c()
cId<-1

#for number of observations, parameters and statistics 
N<-c()
G<-c()
Beta0 <- c()
Beta1 <- c()
Beta2 <- c()
Icc<-c()
Mor<-c()
Soi<-c()
Ior_80lower<-c()
Ior_80upper<-c()
rmse_ICC<-c()
rmse_MOR<-c()
rmse_SOI<-c()
IOR_80R.OVERLAP<-c()
SSIZE<-c()

#beginning of simulations
# note the opening curly bracket

{
  #prompting user to enter values
  number_of_clusters <- readline(prompt = "Enter number of clusters: ");
  cluster_size <- readline(prompt = "Enter cluster size: ");
  l2_variance<-readline(prompt = "Enter the assumed fixed level 2 variance: ");
  
  actial.ior.lower <-readline(prompt = "Enter the value of actial.ior.lower: ");  
  actial.ior.upper<-readline(prompt = "Enter value of actial.ior.upper:");
  actual.icc<-readline(prompt = "Enter the value actual.icc: ");
  actual.mor<-readline(prompt = "Enter the value of actual.mor: ");
  actual.soi<-readline(prompt = "Enter the value of the actual.soi: ");
  
  reps<-readline(prompt = "How many times do you want to replicate?: ")
  
  
  
  
  cluster_size<-as.integer(cluster_size)
  number_of_clusters<-as.integer(number_of_clusters)
  l2_variance<-as.numeric(l2_variance)
  reps<-as.integer(reps)
  
  
  #filenames data
  data_index_name<-"data_cluster"
  mean_index_name<-"mean_cluster"
  fitted_index_name<-"fitted_cluster"
  size_cluster_name<-paste("size",cluster_size, collapse = " ", sep = "_")
  var<-paste("var",l2_variance, collapse = " ", sep = "_")
  ext<-".csv"
  
  #file names
  data_filename<-paste(data_index_name,number_of_clusters, size_cluster_name,var, ext,
                       collapse = " ", sep = "_")
  
  mean_filename<-paste(mean_index_name,number_of_clusters, 
                       size_cluster_name,var, ext, collapse = " ", sep = "_")
  
  fitted_results_filename<-paste(fitted_index_name,number_of_clusters, size_cluster_name,var,
                                 ext, collapse = " ", sep = "_")
  
  
  for (i in 1:reps) {
    
    
    for (clust_ID in 1:number_of_clusters) {
      u<-rnorm(1, 0,sqrt(l2_variance))
      for (j in 1:cluster_size) {
        x1<-rnorm(1,0,1)
        x2<-rnorm(1,0,1)
        pi<-0.5-1.5*x1 + 0.3*x2+ u
        prob <- 1/(1+exp(-(pi)))
        y <- rbinom(1,1, prob)
        X1<-append(X1,x1)
        X2<-append(X2,x2)
        U<-append(U,u)
        Y<-append(Y,y)
        clusterId<-append(clusterId,cId)
        if(j%%cluster_size==0){
          cId<-cId + 1
        }
      }
    }
    #generated data dataframe
    dat<-data.frame(clusterId,Y,X1,X2,U)
    
    #writing generated data to csv
    write.csv(dat, data_filename)
    #write.csv(dat,"data_cluster_100_size_5_var_1.0.csv")
    
    # print(dat)
    l1size<-j
    l2size<-clust_ID
    ssize<-l1size*l2size
    
    #variance<-0.5
    #residual_var = 3.29
    #actial.ior.lower <-0.046
    #actial.ior.upper<- 0.595
    #actual.icc<-0.132
    #actual.mor<-1.958
    #actual.soi<-0.036
    
    
    variance<-l2_variance
    residual_var = 3.29
    
    actial.ior.lower<-as.numeric(actial.ior.lower)
    
    actial.ior.upper<-as.numeric(actial.ior.upper)
    
    actual.icc<-as.numeric(actual.icc)
    
    actual.mor<-as.numeric(actual.mor)
    
    actual.soi<-as.numeric(actual.soi)
    
    #for (i in 1:length(clusterId)){
    m2<-glmer(Y~1+X1+X2+(1|clusterId),data=dat, family =  binomial (link = "logit"))
    fixed<-summary(m2)
    #print(fixed)
    beta0<- fixed$coefficients[1,1]
    beta1<- fixed$coefficients[2,1]
    beta2<- fixed$coefficients[3,1]
    #beta0<-beta0.append(beta0)
    icc<-function(dat){
      b<-VarCorr(dat)
      residual_var = 3.29# attr(b,'sc')^2
      intercept_var = attr(b$clusterId,'stddev')[1]^2
      return(icc = intercept_var/(intercept_var+residual_var))
    } #computes the ICC
    ICC<-icc(m2)
    #print(ICC)
    rmse_icc<-rmse(actual.icc,ICC)
    
    #MOR
    
    mor<-function(dat){b<-VarCorr(dat)
    residual_var = 3.29# attr(b,'sc')^2
    intercept_var = attr(b$clusterId,'stddev')[1]^2
    return(mor = exp(sqrt(2*intercept_var)*qnorm(0.75)))
    }# computes the MOR
    MOR<-mor(m2)
    rmse_mor<-rmse(actual.mor,MOR)
    
    
    #SOI
    
    soi<-function(dat){b<-VarCorr(dat)
    residual_var = 3.29# attr(b,'sc')^2
    intercept_var = attr(b$clusterId,'stddev')[1]^2
    beta = fixef(dat)[2]
    beta1 = fixef(dat)[3]
    return(soi = pnorm((beta-beta1) / sqrt(2*intercept_var)))
    }# computes the SOI
    SOI<-soi(m2)
    rmse_soi<-rmse(actual.soi,SOI)
    
    
    #Interval Odds Ratio (IOR-80)
    
    ior_80lower<-function(dat){b<-VarCorr(dat)
    residual_var = 3.29# attr(b,'sc')^2
    intercept_var = attr(b$clusterId,'stdd')[1]^2
    beta = fixef(dat)[2]
    beta1 = fixef(dat)[3]
    return(ior_80lower=exp((beta-beta1) +
                             sqrt(2*intercept_var)*qnorm(0.10)))
    }# computes the IOR_80Lower
    IOR_80l<-ior_80lower(m2)
    rmse_ior_l<-rmse(actial.ior.lower,IOR_80l)
    
    
    ior_80upper<-function(dat){b<-VarCorr(dat)
    intercept_var = attr(b$clusterId,'stdd')[1]^2
    residual_var = 3.29# attr(b,'sc')^2
    beta = fixef(dat)[2]
    beta1 = fixef(dat)[3]
    return(ior_80upper=exp((beta-beta1)+
                             sqrt(2*intercept_var)*qnorm(0.90)))
    }# computes the IOR_80Upper
    IOR_80u<-ior_80upper(m2)
    rmse_ior_u<-rmse(actial.ior.upper,IOR_80u)
    
    # Computing %ge Overlap in IOR-80
    
    if(actial.ior.upper >= IOR_80l)
    {
      IOR80R.overlap <- T
      T1 <-IOR80R.overlap
    }else{
      IOR80R.overlap <- F
      T1 <-IOR80R.overlap
    }
    if(actial.ior.lower <IOR_80u && actial.ior.lower >=IOR_80l)
    {
      IOR80R.overlap <- T
      T2<-IOR80R.overlap
    }else{
      IOR80R.overlap <- F
      T2<-IOR80R.overlap
    }
    if (T1== T && T2==T)
    { IOR80R.overlap<- T
    }else{
      IOR80R.overlap<-F}
    
    #Fitted values
    Beta0<-append(Beta0,beta0)
    Beta1<-append(Beta1,beta1)
    Beta2<-append(Beta2,beta2)
    Icc<-append(Icc,ICC)
    Mor<-append(Mor,MOR)
    Soi<-append(Soi,SOI)
    Ior_80lower<-append(Ior_80lower, IOR_80l)
    Ior_80upper<-append(Ior_80upper, IOR_80u)
    rmse_ICC<-append(rmse_ICC, rmse_icc)
    rmse_MOR<-append(rmse_MOR, rmse_mor)
    rmse_SOI<-append(rmse_SOI, rmse_soi)
    IOR_80R.OVERLAP<-append(IOR_80R.OVERLAP, IOR80R.overlap)
    SSIZE<-append(SSIZE, ssize)
    N<-append(N, l1size)
    G<-append(G, l2size)
  }
  # fitted values's dataframe
  calc_results<-data.frame(Beta0, Beta1, Beta2, Icc, Mor, Soi, 
                           Ior_80lower, Ior_80upper, rmse_ICC, 
                           rmse_MOR, rmse_SOI, IOR_80R.OVERLAP, SSIZE,
                           N, G)
  # saving fitted values to csv
  write.csv(calc_results, fitted_results_filename)
 
  
  #fitted mean values--rounde to 4 decimal places
  mean_beta0<-round(mean(Beta0),4)
  mean_beta1<-round(mean(Beta1),4)
  mean_beta2<-round(mean(Beta2),4)
  mean_Icc<-round(mean(Icc),4)
  mean_Mor<-round(mean(Mor),4)
  mean_Soi<-round(mean(Soi),4)
  mean_Ior_80lower<-round(mean(Ior_80lower),4)
  mean_Ior_80upper<-round(mean(Ior_80upper),4)
  mean_rmse_ICC<-round(mean(rmse_ICC),4)
  mean_rmse_MOR<-round(mean(rmse_MOR),4)
  mean_rmse_SOI<-round(mean(rmse_SOI),4)
  
  #mean fitted values's dataframe
  mean_results<-data.frame(mean_beta0, mean_beta1, mean_beta2, mean_Icc, mean_Mor,
                           mean_Soi, mean_Ior_80lower, mean_Ior_80upper, mean_rmse_ICC, 
                           mean_rmse_MOR, mean_rmse_SOI)
  
  #fitted mean results saved to csv
  write.csv(mean_results, mean_filename)
}
