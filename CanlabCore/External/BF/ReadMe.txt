-----------------------------------------------------------------------------------

BAYES FACTORS - Last updated: 08/12/2014
-----------------------------------------------------------------------------------

Matlab functions for calculating default Bayes Factors (BF10) quantifying the evidence in favour of H1 over H0. Credit for the intellectual development of these algorithms goes to the original authors. I simply implemented them in Matlab. Special thanks go to EJ Wagenmakers for helping fix a programming issue.

Sam Schwarzkopf, UCL
-----------------------------------------------------------------------------------

corrbf(r,n)

BF10 for a correlation with Pearson's r and sample size n.
(See Wetzels et al., 2012 for details)
-----------------------------------------------------------------------------------

corrbfrep(r,n)

BF10 for replicating a correlation with Pearson's r and sample size n.
Uses a uniform prior between 0-1 to test if correlation goes in same direction.
If original r was negative, multiply r with -1.
(See Boekel et al., 2014 for details)
-----------------------------------------------------------------------------------

t1smpbf(t,n,r)

BF10 for a one-sample t-test with t-statistic t and sample size n. The additional input r is the scale factor. This defaults to 0.707 as is the case at Rouder's website (http://pcl.missouri.edu/bf-one-sample). 
(See Rouder et al., 2009 for details)
-----------------------------------------------------------------------------------

t2smpbf(t,nx,ny,r)

BF10 for a two-sample t-test with t-statistic t and sample sizes nx and ny. The additional input r is the scale factor. This defaults to 0.707 as is the case at Rouder's website (http://pcl.missouri.edu/bf-two-sample). 
(See Rouder et al., 2009 for details)
-----------------------------------------------------------------------------------

binombf(k,n,p)

BF10 for a binomial test with k successes, n trials overall, and base probability p.   This is based on the Bayes Factor example on Wikipedia (http://en.wikipedia.org/wiki/Bayes_factor). It should be the same as the default prior assumed by Rouder's web applet (http://pcl.missouri.edu/bf-binomial).

 
