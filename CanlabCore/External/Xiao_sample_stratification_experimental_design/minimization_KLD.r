###########################################################################################
# Date	:	12/10/2010						                              #
# Author	:	LX							                              #
# Program	:	randomization_KLD.R					                        #
# Description:	Calculates the KLD index D based on Endo et al. (2006)	            #
# Procedure for allocation:  1). The first 2*ntrt subjects are randomly allocated to each #
#                                group with two subjects per group.                       #
#                            2). if the number of difference in the assigned sample sizes #
#                                between two groups >=Dn, then new subject is allocated to#
#                                the group with least subjects with probability p_Dn.     #
#                            3). Compares the D for treatments and control and assigns    #                          
#                                descendingly Pks  to the groups with ascendingly ordered # 
#                                imbalance scores and equal probability to the groups     #   
#                                with equal imbalance score.                              #                   
#                                                                                         # 
###########################################################################################


##################################################################################################
# Input data: Three separated CVS formatted files under the same folder: Continuous variables 
#             (continuous_covariates) , categorical variables (categorical_covariates), and
#             group assigned ( group_old) if there are subjects already randomized whenever
#             nassigned^=0. The input parameters include the number of treatments (ntrt), 
#             the number of assigned subjects (nassigned), the file name of output allocations
#             (outfile), the maximum allowable difference in total number of subjects 
#             between any two groups (Dn) and probability p_Dn, and  a vector of probability for ascendingly 
#             ordered imbalance scores(Pk).  
#               
# Note: 1). All subjects which include the assigned and new subjects need to be included in 
#       the files (continuous_covariates and categorical_covariates). 
#       2). each file includes unique subject identifier (i.e. studyid) as first column.                          
###################################################################################################
# Example parameter input:
# folder	=	"c:\\randomization\\"
# ntrt = 3
# continuous covariates	=	"continuous.csv"
# categorical covariates=       "categorical.csv"
# nassigned            =        5
# group_old                 =   "group_assigned.csv"
# outfile	=	"group_new.csv"
# Dn = 4
# p_Dn =0.9
# Pk = c(2/3,1/6,1/6) 
# Note: "group_assigned.csv" as input for group_old is the outfile from previous assignment and 
#        should include 5 assigned subjects in this example due to nassigned=5.
#       "group_new.csv" as an outfile will include already assigned and new assignments.  
###################################################################################################

# Yoni's notes:  Pk can be .9 to favor balance while still allowing some randomness
# Efron [13] proposed Dn = 1 and p Dn = 2/3.  XL recommends pDn of .8 or .9


random.allocation	<-	function(folder,ntrt, continuous_covariates, categorical_covariates, nassigned, group_old, outfile, Dn, p_Dn, Pk) {
  dir			<-  	folder
  confile		<-	paste(dir,continuous_covariates,sep="/")
  con<-read.csv(confile,header=T, sep=',')
  
  
  catfile             <-      paste(dir,categorical_covariates,sep="/")
  cat<-read.csv(catfile,header=T, sep=',')
  
  if (nassigned!=0){
    groupfile             <-      paste(dir,group_old,sep="/")
    group_old<-read.csv(groupfile,header=T, sep=',')}
  
  n<-nrow(con)
  n_new<-n-nassigned
  
  if (nassigned == 0) {
    group <- rep(-1, n_new)
    random_number <-
      rep(-1, n_new)
  } else{
    group <- c(group_old[, 2], rep(-1, n_new))
    random_number <- c(group_old[, 3], rep(-1, n_new))
  } 
  
  
  ###########################################################################################################
  # First 2* ntrt subjects are randomly allocated the groups at two subjects per group.
  ###########################################################################################################
  
  
  if (nassigned==0) {
    seed<-rnorm(1,0,10000)
    random_number[nassigned+1]<-seed
    set.seed(seed)
    permute<-sample(c(1:ntrt,1:ntrt),2*ntrt,replace=F)
    group[nassigned+1]<-permute[1]
    nassigned <- nassigned + 1}
  while (nassigned >= 1 &
         nassigned < (2 * ntrt) & n >= nassigned) {
    seed <- random_number[1]
    random_number[nassigned + 1] <- seed
    set.seed(seed)
    permute <- sample(c(1:ntrt, 1:ntrt), 2 * ntrt, replace = F)
    group[nassigned + 1] <- permute[nassigned + 1]
    nassigned <- nassigned + 1
  }
  
  #######################################################################################################################################################
  #(2*ntrt+1)th subject and after, assign to the treatment arm with fewest number of subject with high probability p_Dn or Pk to the ascending KLD index.
  #######################################################################################################################################################
  
  nc<-ncol(con)-1
  d0_con<-matrix(0,nc,ntrt)
  c<-nassigned+1
  groupc<-group[1:c]
  conc<-con[1:c, ]
  colc<-cat[1:c, ]
  d0<-numeric(ntrt)
  
  while (c>(2*ntrt) & c<=n)  {
    len<-numeric(ntrt) 
    
    # compute N of each tx group
    for (l in 1:ntrt) {
      len[l]<-length(group[1:(c-1)][group[1:(c-1)]==l])
    }
    
    # rand gen
    seed<-rnorm(1,0,10000)
    random_number[c]<-seed
    set.seed(seed)
    
    # if group N imbalance exceeds Dn
    if ((max(len)-min(len))>=Dn) { 
      p<-(1-p_Dn)/(ntrt-1)
      select<-sample(sort(len),1,prob=c(p_Dn,rep(p,(ntrt-1)))) # select from smaller-N group w/ p_Dn, larger groups w/ 1-p_Dn. will select the N..?
      
      if (unique(select)==min(len)){ # if we "passed" p_Dn threshold, and should assign S to smaller N group
        
        # allocate Subject.  don't understand exactly
        select1<-(1:length(len))[len==select]
        if (length(select1)==1) {
          group[c]<-select1
        } 
        else {
          group[c]<-sample(select1,1)
        }
      }
    }
    
    # subject has not been allocated. compute KLDs for all groups, save into 'd'
    if (group[c]==-1)    {
      
      # continuous
      for (i in 2:ncol(con)) {
        for (a in 1:ntrt){
          for (t in 1:ntrt){
            groupc[c]<-a 
            if (t==a){temp<-0} else {    
              # if no var in continuous cov, throws error
              temp<- -2 + 0.5*( (mean(conc[ ,i][groupc==a])-mean(conc[ ,i][groupc==t]))**2 + (var(conc[ ,i][groupc==a])+var(conc[ ,i][groupc==t])))  *(1/var(conc[ ,i][groupc==a]) + 1/var(conc[ ,i][groupc==t]) )                         
              d0_con[i-1,a]<-d0_con[i-1,a]+temp }}  
        }
      }
      
      # ?
      nca<-ncol(cat)-1
      d0_col<-matrix(0,nca,ntrt)
      
      # categorical
      for (i in 2:ncol(cat)) { 
        for (a in 1:ntrt){
          
          for (t in 1:ntrt){  
            for (j in (unique(cat[ ,i]))){
              groupc[c]<-a    
              if (t==a){temp0<-0} else {       
                temp0 <- (length(colc[ ,i][groupc==a & colc[ ,i]==j ])/length(colc[ ,i][groupc==a])-length(colc[ ,i][groupc==t & colc[ ,i]==j ])/length(colc[ ,i][groupc==t]))*(log(length(colc[ ,i][groupc==a & colc[ ,i]==j ])/length(colc[ ,i][groupc==a])+1e-12)-log(length(colc[ ,i][groupc==t & colc[ ,i]==j ])/length(colc[ ,i][groupc==t])+1e-12))
                d0_col[i-1,a]<-d0_col[i-1,a]+temp0
              }
            }
          }
        }
      }
      
      # compute final d (cov imbalance) measures for all groups. in d0, columns are groups, rows are covs
      for (t in 1:ntrt){
        d0[t]<-sum(d0_con[ ,t])+sum(d0_col[ ,t])
      }
      d<-d0
      
      # select group using Pk probabilities over d imbalances per group
      seed<-rnorm(1,0,10000)
      random_number[c]<-seed
      set.seed(seed)
      select<-sample(sort(d),1, prob=Pk)
      select1<-(1:length(d))[d==select]
      if (length(select1)==1) {group[c]<-select1} else {
        group[c]<-sample(select1,1)} # the allocation
    }
    
    # increment values
    nassigned<-nassigned+1
    c<-nassigned+1
    groupc<-group[1:c]
    conc<-con[1:c, ]
    colc<-cat[1:c, ]
    
  }
  
  
  studyid<-con[ ,1] # grab studyid from convenient place
  out<-data.frame(studyid,group,random_number)
  
  write.csv(out,paste(dir,outfile,sep="/"),row.names=F,)
}
