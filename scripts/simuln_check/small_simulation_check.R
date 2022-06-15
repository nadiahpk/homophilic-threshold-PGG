rm(list=ls())

# alternative to standard sample() function that always treats x as a vector to
#  sample from
resample <- function(x, ...) x[sample.int(length(x), ...)]

# model parameters (these are for Fig. 3b)
n = 8
tau = 4
W = 2
X = -0.5
Y = 3
Z = 0

# range of homophily values (h) to consider
n_h = 11
hs = seq(0,1,length.out=n_h)
# range of initial co-operator frequencies (p0) to consider
n_p0 = 11
p0s = seq(0.01,0.99,length.out=n_p0)

# matrix to store all the results for final co-operator frequency (p)
all_ps = matrix(NA,n_h,n_p0)

# allowed models are currently 'leader-driven' and 'members-recruit'
model = 'members-recruit'

time0 = proc.time()

# loop over all homophily values and all initial co-operator frequencies
for ( i in 1:n_h )
  {
  h = hs[i]
  for ( j in 1:n_p0 )
    {
    # initialise the simulation
    p0 = p0s[j]
    cat('**** h = ',h,', p0 = ',p0,'\n',sep='')
    p = p0
    dt = 0.01
    # t_max needs to be high enough to get away from transient dynamics, but
    #  not too high because this is a stochastic model and eventually we'll hit
    #  one of the absorbing states (p=0 or p=1), which is not what we're looking
    #  for
    t_max = 20
    # higher values of n_rep should give more accurate results--more groups
    #  will be sampled in each timestep
    n_rep = 100
    t = 0
    ts = seq(0,t_max,dt)
    ps = rep(NA,length(ts)+1)
    ps[1] = p
    
    # loop over time
    for ( t_i in 1:length(ts) )
      {
      t = ts[t_i]
      
      # these variables are used to calculate the average co-operator and
      # defector fitnesses
      tot_fit_C = 0
      tot_n_C = 0
      tot_fit_D = 0
      tot_n_D = 0

      # loop over number of groups to be sampled
      for ( k in 1:n_rep )
        {
        # initialise the group with just the first individual
        ids = c(1,rep(NA,n-1))
        next_id = 2
        # now select the rest of the individuals
        for ( l in 2:n )
          {
          r = runif(1)
          
          if ( model == 'leader-driven' )
            {
            if ( r < h )
              {
              # recruit kin
              ids[l] = 1
              } else
              {
              # recruit non-kin
              ids[l] = next_id
              next_id = next_id + 1
              }
            } else if ( model == 'members-recruit' )
            {          
            if ( r < h )
              {
              # recruit kin
              ids[l] = resample(ids[1:(l-1)],1)
              } else
              {
              # recruit non-kin
              ids[l] = next_id
              next_id = next_id + 1
              }
            }
          }

        # choose strategies for each family
        #  (next_id-1 is the number of families)
        strategies = sample(c('C','D'),next_id-1,replace=T,prob=c(p,1-p))
        # compute number of co-operators
        n_C = sum(strategies[ids]=='C')

        # assign fitnesses to each individual in the group depending on
        #  strategies and whether the threshold was met
        threshold_met = (n_C>=tau)
        if ( threshold_met )
          {
          tot_fit_C = tot_fit_C + W*n_C
          tot_n_C = tot_n_C + n_C
          tot_fit_D = tot_fit_D + Y*(n-n_C)
          tot_n_D = tot_n_D + n - n_C
          } else
          {
          tot_fit_C = tot_fit_C + X*n_C
          tot_n_C = tot_n_C + n_C
          tot_fit_D = tot_fit_D + Z*(n-n_C)
          tot_n_D = tot_n_D + n - n_C
          }
        }

      # calculate average fitnesses of co-operators and defectors in this
      #  timestep, based on all the groups we sampled above
      
      # if we didn't sample any co-operators, then we can assume they'll almost
      #  always be playing against defectors
      fit_C = ifelse(tot_n_C==0,X,tot_fit_C/tot_n_C)
      # if we didn't sample any defectors, then we can assume they'll almost
      #  always be playing against co-operators
      fit_D = ifelse(tot_n_D==0,Y,tot_fit_D/tot_n_D)
      
      # do the adaptive dynamics
      p = max(0,min(1,p + p*(1-p)*(fit_C-fit_D)*dt))
      
      # store the results for this timestep
      ps[t_i+1] = p
      }
    cat('p_final = ',p,'\n',sep='')
    # store the final co-operator frequency for this simulation to the matrix
    all_ps[i,j] = p
    }
  }

time1 = proc.time()
print(time1-time0)

# draw the bifurcation diagram
ppp = rep(as.vector(p0s),each=n_p0)
iii = which(ppp>0 & ppp<1)
plot(rep(as.vector(hs),n_p0)[iii],as.vector(all_ps)[iii],xlab='h',ylab='p',pch=19,col='red')

# write to csv
df <- data.frame(cbind(hs, all_ps))
colnames(df) <- c('h', p0s)
write.csv(df,'../../results/simuln_check/small_simulation_check.csv', row.names = FALSE)
