import numpy as np

# add a new family with abundance 1
def add_new(abunV, famsV):

    # add new family
    if 1 in abunV:

        idx = abunV.index(1)
        famsV[idx] += 1

    else:

        abunV.append(1)
        famsV.append(1)

    return(abunV, famsV)

# recruit a new member to existing family
def recruit_family(abunV, famsV):

    # create the cumulative sum
    cumuV = [ abunV[0]*famsV[0] ]
    for j in range(1, len(abunV)):

        cumuV.append( cumuV[j-1] + abunV[j]*famsV[j] )

    # choose a recruitor
    rec_idx = np.random.randint(1, cumuV[-1]+1) # index in [1, group size]

    # find the abundance of recruitor's family, and the no. of families with that abundance

    cumu_idx = 0
    while rec_idx > cumuV[cumu_idx]:
        cumu_idx += 1

    rec_abun = abunV[cumu_idx] # abundance of recruitor's family
    rec_fams = famsV[cumu_idx] # number of families with that abundance

    # one of those families now has +1 individuals
    new_abun = rec_abun + 1

    # add new-abundance family

    if new_abun in abunV:

        # if families with that abundance already exist in the gorup, 
        # add one to the count of such families

        idx = abunV.index(new_abun)
        famsV[idx] += 1

    else:

        # if families with that abundance don't yet exist in the group, add it

        abunV.append(new_abun)
        famsV.append(1)

    # remove recruitor-abundance family 

    if rec_fams == 1:

        # if the recruitor's family was the only family with that abundance,
        # we need to remove this abundance category altogether

        del abunV[cumu_idx]
        del famsV[cumu_idx]

    else:

        # if there are other families with the same abundance as recruitor's,
        # reduce the count of such families by 1

        famsV[cumu_idx] -= 1

    return(abunV, famsV)

