# allocate potentially eligible Ps to study.
# steps:  copy info for new patients to study_allocated.csv, run this script, then re-open file to see their allocation
# 1 = Placebo, 2 = PRT

# initial set up.
folder = '/Users/yoni/Dropbox/Placebo injection/Scripts and Analyses/allocation'
setwd(folder)
source("minimization_KLD.R")

n_trt = 2
n_assigned = 0
myDn = 2
myp_Dn = .5 # in her and my simulations, .9 and 1 (or i am doing .95) for the two vals below yielded best results
myPk = c(0.95, 0.05)

# read in already allocated Ps. 
alloc_file = 'study_allocated.csv'
alloc = read.csv(alloc_file, sep=',')

# save to back up folder just in case, b/c will be overwritten
write.csv(alloc, file = paste('bak/', 'study_allocated', Sys.time(), '.csv'), row.names=FALSE)

# separate into two cov files and save. 
covs = read.csv('study_covs.csv', sep=',')
write.csv(covs[, c(1, 2, 4)], file = 'tmp/con_cov.csv', row.names=FALSE)
write.csv(covs[, c(1, 3, 5)], file = 'tmp/cat_cov.csv', row.names=FALSE)

# allocate new Ss, saving to same file (overwriting)
random.allocation(folder = folder, ntrt= n_trt, continuous_covariates='tmp/con_cov.csv', categorical_covariates='tmp/cat_cov.csv', nassigned=nrow(alloc), group_old = alloc_file, outfile = alloc_file, Dn=myDn, p_Dn=myp_Dn, Pk= myPk)

# check current allocation values
alloc = read.csv(alloc_file, sep=',') # re-read w/ new values
sum(alloc$group==1)
sum(alloc$group==2)

# check covariate balance
dat = merge(alloc, covs)
boxplot(age ~ group, data=dat)
title('age')
boxplot(pain ~ group, data=dat)
title('baseline pain')

t.test(dat$pain ~ dat$group, paired = FALSE)
t.test(dat$age ~ dat$group, paired = FALSE)

tabl=table(dat$gender, dat$group)
tabl
chisq.test((tabl))

tabl=table(dat$opioid, dat$group)
tabl
chisq.test((tabl))

options(max.print=1000000)
alloc
