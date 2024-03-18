r"""
Code dealing with conjecture 1.9

Summary goes here

Copyright stuff is here because I hosted it on github as well, and they insist on something like that

EXAMPLES::

<Lots and lots of examples>

AUTHORS:

- Tucker Ervin (2020-10-13): initial version

- Blake Jackson (2020-10-13): initial version

"""

# ****************************************************************************
#       Copyright (C) 2020 Tucker Ervin & Blake Jackson tjervin@crimson.ua.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

def rCancel(rSeq):
    r'''
    This function takes in a sequence of indices and removes duplicates next to each other

    Input:

    - ''rSeq''  -- sequence of indices; Represents the indices of the product ''r_{i_1}r_{i_2}\\dots r_{i_n}''

    Output: ''red'' -- reduced sequence of indices
    '''
    red = rSeq
    j = 0 # Begins at the first element of the sequence
    while j < (len(red)-1):
        if red[j] == red[j+1]: # Checks if the current element is equal to the next element
            red1 = red[0:j]
            red2 = []
            if (j+2 <= len(red)):
                red2.extend(red[j+2:])
            red = red1 + red2 # Removes the two equal neighboring elements
            j = 0 # Begins again
        else:
            j += 1 # Advances one element
    return red


def mMutation(M, w):
    r'''
    This function takes in a matrix and outputs its image under matrix mutation

    Input:
    - ''M'' -- matrix;

    - ''w'' -- sequence of mutations; The mutations range from 1 to rank(M).

    Output: ''mut'' -- mutated Matrix
    '''
    length = len(w)
    r = M.nrows()
    c = M.ncols()
    mut = copy(M) # Builds a a copy of the original matrix
    if length == 0: # In case no mutation happens
        return mut
    else:
        for i in w:
            if i < 1 or i > r or i > c: # Checks to see if we are mutating at non-vertices
                print("Invalid mutation at vertex: ", i)
                return mut
    if length == 1: # Only one mutation happens
         k = w[0]
         for i in range(r):
             for j in range(c): # Standard matrix mutation here
                 if i == k-1 or j == k-1: # Has the -1's here to account for sage indexing at 0
                     mut[i,j] = -M[i,j]
                 elif M[i,k-1]*M[k-1,j] > 0: # Has the -1's here to account for sage indexing at 0
                     mut[i,j] = M[i,j] + M[i,k-1]*M[k-1,j]*(M[i,k-1]/abs(M[i,k-1])) # Has the -1's here to account for sage indexing at 0
                 else:
                     mut[i,j] = M[i,j]
    else: # Multiple mutations happen
        mut = mMutation(mMutation(M,w[0:length-1]), [w[length-1]]) # Recursively "goes up" the mutation sequence
    return mut

def corrGIM(M, l):
    r'''
    This function takes in a matrix and a sequence representing a total order, returning a GIM

    Input:
    - ''M'' -- skew-symmetrizable matrix;

    - ''l'' -- sequence representing linear order; the length of the linear order equals the rank of M

    Output: ''A'' -- appropriate GIM
    '''
    n = M.nrows()
    length = len(l)
    A = zero_matrix(n,n) # Starts A as a zero matrix
    if n != length:
        print("Invalid order")
    else:
        for i in range(n):
            for j in range(n): # Standard rules for constructing GIM
                if i == j:
                    A[i,j] = 2
                else:
                    li = l.index(i+1) # Finds the index of the first and only instance of i+1 in the order
                    lj = l.index(j+1) # Ditto ^
                    if li < lj: # Means that i+1 < j+1 in the order
                        A[i,j] = M[i,j]
                    else:
                        A[i,j] = -M[i,j]
    return A

def rMutation(i, bM, cM, w):
    r'''
    This function takes in an index i, the initial exchange matrix, the initial c-vector matrix, and a mutation sequence w
    Reduces in a way which does not well define the g_i^w

    Input:
    - ''i'' -- index; Indexing begins at 1 not 0, as we are accustomed to

    - ''bM'' -- initial exchange matrix; Both matrices should be square

    - ''cM'' -- initial c-vector matrix; Both matrices should be square

    - ''w'' -- mutation sequence;

    Output: ''rW'' -- an index sequence representing ''r_i^w''
    '''
    wPrime = []
    length = len(w)-1
    rW = []
    if i < 1 or i > bM.nrows(): # Checks to see if this is a valid index we are picking
         print ("Invalid mutation at index: ", i)
         return w
    if length == -1: # When w is the empty mutation
         rW = [i]
    elif length == 0: # This is when w is a single mutation
         k = w[0]
         if (i == k): # Identity we showed
             rW = [i]
         else:
             cK = bM[i-1,k-1] * cM[k-1, :] # Has the -1's here to account for sage indexing at 0
             if cK > 0: # Definition of r_i^[k] from paper
                 rW = [k, i, k]
             else:
                 rW = [i]
    else:
         wPrime = w[0:length]
         k = w[length] # Hence w = wPrime[k]
         if (i == k): # Identity we showed
             rW.extend(rMutation(i, bM, cM, wPrime)) # Recursive step here
         else:
             M = mMutation(block_matrix([[bM,cM]]), wPrime) # Puts the matrix in [ B C ] form so we can perform matrix mutation on it
             # The cluster algebra mutation provided with sage does column c-vectors instead, and it was messing up my calculations for some reason
             r = M.nrows()
             cK = M[i-1,k-1] * M[k-1, r:] # Has the -1's here to account for sage indexing at 0
             rW1 = rMutation(i,bM,cM,wPrime)
             rW2 = [] # Corresponds to the case when r_i^(w'[k]) = r_i^(w')
             if cK > 0:
                       rW2.extend(rMutation(k,bM,cM,wPrime)) # Corresponds to the case when r_i^(w'[k]) = r_k^(w')r_i^(w')r_k^(w')
             rW = rW2 + rW1 + rW2 # Cocatenates the lists
    rW = rCancel(rW) # Automatically removes some of the r_i's when possible
    return rW


def rAction(iSeq, aM):
    r'''
    This function takes in an index sequence representing r_i^w and returns a matrix representing its action under pi

    Input:
    - ''iSeq'' -- index sequence; Indices begin at 1 not 0

    - ''aM'' -- GIM from linear ordering;

    Output: ''A'' -- square matrix representing '\\pi(r_i^w) acting on elements of '\\Gamma' on the left
    '''
    r = aM.nrows()
    c = aM.ncols()
    A = matrix(r,c,0) # Gives a fresh matrix to represent the action
    length = len(iSeq)
    if length == 0:
        return A
    i = iSeq[length-1] # As the r_i's act on the left, we need to evalute the last one first
    if length == 1:
        for j in range(r):
            for k in range(c): # This follows the rules set for the action of r_i given in paper
                if j == k:
                    if k == i-1: #-1 here to account for indexing
                        A[j,k] = -1
                    else:
                        A[j,k] = 1
                elif j == i-1:
                    A[j,k] = -aM[k,j]
                else:
                    A[j,k] = 0
    else:
        A = rAction(iSeq[0:length-1], aM)*rAction([i], aM) # Recursively multiplies on the left
    return A


def rEqual(iSeq, jSeq, aM):
    r'''
    This function checks if the actions of two sequences representing 'r_i^w' and 'r_j^v' are equal

    Input:
    - ''iSeq'' -- index sequence; Indexing begins at 1 not 0

    - ''jSeq'' -- index sequence;

    - ''aM'' -- GIM from linear ordering;

    Output: ''result'' -- boolean
    '''
    A1 = rAction(iSeq, aM) # Gets the two actions to compare
    A2 = rAction(jSeq, aM)
    result = (A1 == A2) # Compares them
    if not result:
        print("Actions are not equal")
        print("At seq. iSeqW")
        print(iSeq)
        print(A1)
        print("At seq. iSeqV")
        print(jSeq)
        print(A2)
    return result

def allREqual(M, aM, w, v):
    r'''
    This function checks if 'r_i^w = r_i^v' for all 'i'

    Input:
    - ''M'' -- 2n x n matrix;

    - ''aM'' -- corresponding GIM;

    - ''w'' -- mutation sequence; Indexing begins at 1 not 0

    - ''v'' -- mutation sequence;

    Output: ''result'' -- boolean
    '''
    n = M.nrows() # Gets the number of indices
    result = true # default return is true
    bM = M[:,:n]
    cM = M[:,n:]
    for i in [1..n]:
        iSeqW = rMutation(i, bM, cM, w) # Calculates the two sequences of r_i's
        iSeqV = rMutation(i, bM, cM, v)
        if not rEqual(iSeqW, iSeqV, aM): # compares the two. If false they are not equal
            print('Linear order is not suitable at index: ', i)
            result = false
            break
    return result

def cEqual(M, w, v):
    r'''
    This function determines if the C-vectors are equal under matrix mMutation

    Input:
    - ''M'' -- n x 3n matrix;

    - ''w'' -- mutation sequence; Indexing begins at 1 not 0

    - ''v'' -- mutation sequence;

    Output: ''result'' -- boolean
    '''
    n = M.nrows()
    cW = mMutation(M, w)[:,n:] # Gets the C matrices
    cV = mMutation(M, v)[:,n:]
    result = (cW == cV) #compares the two
    return result

def conjSatisfied(M, B, w, v, l):
    r'''
    This function checks to see if the conjecture is satisfied for two given mutation sequences

    Input:
    - ''M'' -- 2n x n matrix;

    - ''B'' -- initial exchange matrix;

    - ''w'' -- mutation sequence; Indexing begins at 1 not 0

    - ''v'' -- mutation sequence;

    - ''l'' -- sequence representing linear order; the length of the linear order equals the rank of bM

    Output: ''result'' -- boolean; Default is true
    '''
    n = M.nrows()
    aM = corrGIM(B, l)
    result = true
    if cEqual(M, w, v):
        result = allREqual(M, aM, w, v)
    else:
        print("C^w does not equal C^v. Does not satisfy assumptions of conjecture")
        result = false
    return result




'''
It's important to remember that we only care about based loops in the exchange graph. Now we need a way to generate loops in the exchange graph. This is probably not an efficent way to do this. adj is an adjacency matrix, i is that starting seed/vertex. Returns a mutation sequence that gives a loop in the exchange graph. Indexing starts at 1 not 0.
'''

def randomWalk(oB):
    if oB.nrows() == oB.ncols():
        n = oB.ncols()
        B = oB.set_immutable() #immutable version of B if we need it
        G = fullExGraph(oB,Graph()) #creates the associhedron for B
        MarkovB = matrix(QQ,G.adjacency_matrix()) #A markov chain for the associhedron
        u = sum(posB.columns()) #this gives a vector of the row sums of our adjacency matrix
        for i in [1..n]:
            for j in [1..n]:
                MarkovB[i-1,j-1] = MarkovB[i-1,j-1]/u[i-1] #turns B into a markov chain matrix with uniform probabilities
        walk = [] #empty walk vector
        walk = rCancel(walk) #this removes any backtracking/extraneous steps
    else:
        print("B matrices must be square.")
    return walk




'''
##############################################################################################
USING B INSTEAD OF M GIVES THE DESIRED EXCHANGE GRAPH.
##############################################################################################
'''

def exGraphMutation(oM, gO):
    r'''
    Adds the vertices and edges connected to a matrix in the labelled exchange graph

    Input:
    - ''oM'' -- 'n \\times 2n';
    - ''gO'' -- graph that is being mutated;

    Output: ''G'' -- Graph with mutations added to it
    '''
    G = copy(gO)
    o = G.order()
    M = copy(oM)
    M.set_immutable()
    n = M.nrows()
    if o == 0:
        G.add_vertex(M)
    else:
        for i in [1..n]:
            mut = mMutation(M, [i])
            mut.set_immutable()
            edge = {M, mut}
            strLabel = '%i'%i
            if G.has_vertex(mut):
                if not G.has_edge(edge):
                    G.add_edge(edge, label=strLabel)
                    #print(str)
            else:
                G.add_vertex(mut)
                G.add_edge(edge, label=strLabel)
                #print(str)
    return G

def exGraph(M, gO):
    r'''
    Adds the vertices and edges connected to a labelled exchange graph for all vertices currently in the graph

    Input:
    - ''M'' -- 'n \\times 2n';
    - ''gO'' -- graph that is being mutated;

    Output: ''G'' -- Graph with mutations added to it
    '''
    G = copy(gO)
    if G.order() == 0:
        G = exGraphMutation(M,G)
    vertices = G.vertices()
    n = M.nrows()
    for vert in vertices:
        deg = G.degree(vert)
        if deg != n:
            newG = exGraphMutation(vert, G)
            G = G.union(newG)
    return G

def fullExGraph(M):
    r'''
    Creates the full labelled exchange graph for a given matrix

    Input:
    - ''M'' -- 'n \\times 2n';

    Output: ''G'' -- Graph with mutations added to it

    WARNING: This function is quite slow. Do not try for anything above 'A_4'
    '''
    G = Graph(0)
    newG = exGraph(M, G)
    while G != newG:
        G = copy(newG)
        newG = exGraph(M, G)
    return G

def getMutationSeq(G, cycle, vert):
    r'''
    Take in a graph and outputs a sequence representing the edges taken in a cycle
    
    Input:
    - ''G'' -- Graph;
    
    - ''cycle'' -- cycle in graph; Assume this is actually a cycle
    
    - ''vert'' -- vertex that the cycle starts at;
    
    Output: ''muSeq'' -- Sequence representing edges taken starting at cycle starting point
    '''
    muSeq = []
    if vert not in cycle:
        print('Vertex was not a member of the cycle')
        return muSeq
    length = len(cycle)
    start = 0
    muStart = []
    muEnd = []
    mu = 1
    for i in [0..length-1]:
        if cycle[i] == vert:
            start = i
            break
        else:
            if i < length-1:
                mu = int(G.edge_label(cycle[i], cycle[i+1]))
            else:
                mu = int(G.edge_label(cycle[i], cycle[0]))
            muEnd.append(mu)
    for i in [start..length-1]:
        if i < length-1:
            mu = int(G.edge_label(cycle[i], cycle[i+1]))
        else:
            mu = int(G.edge_label(cycle[i], cycle[0]))
        muStart.append(mu)
    muSeq = muStart + muEnd
    return muSeq

def getMinSeq(startVert, G, endVert):
    minPath = []
    if startVert == endVert:
        return minPath
    rG = G.random_spanning_tree(output_as_graph=True)
    paths = rG.all_paths(startVert, endVert)
    minPath = paths[0]
    minMut = []
    length = len(minPath)
    for i in [0..length-2]:
        mu = 1
        mu = int(G.edge_label(minPath[i],minPath[i+1]))
        minMut.append(mu)
    return minMut

def conjSatisCycle(M, G, l):
    r'''
    Take in a LEG and outputs whether or not the conjecture is satisfied for cycles containing M
    
    Input:
    - ''M'' -- Matrix that starts the EGG;
    
    - ''G'' -- Graph;
    
    - ''l'' -- linear order;
    
    Output: ''result'' -- boolean;
    '''
    result = false
    n = M.nrows()
    B = M[:, :n]
    
    
    return result

def lVector(r, n, A):
    r'''
    Calculates the L vectors given a r_i^w
    
    Input:
    - ''r'' -- sequence representing r_i^w;
    - ''n'' -- rank of the cluster algebra;
    - ''A'' -- GIM;
    
    Output: ''l'' -- vector corresponding to r_i^w:
    '''
    m = len(r)
    index = floor(m/2)
    i = r[index]
    a = [0 for j in [1..n]]
    a[i-1] = 1
    a = vector(a)
    if index > 0:
        #print("r_" + str(i) + "^w is " + str(r))
        #print("g_" + str(i) + "^w is " + str(r[0:index]))
        M = rAction(r[0:index],A)
    else:
        #print("r_" + str(i) + "^w is " + str(r))
        M = identity_matrix(n)
        
    #if M != identity_matrix(n):
        #print("Action of g_" + str(i) + "^w: ")
        #print(M)
        #print("Action of s_" + str(i) + ": ")
        #print(rAction([i],A))
    l = M * a
    
    return l

def lMatrix(r, A):
    r'''
    Calculates the L^w matrix given a list of  r_i^w
    
    Input:
    - ''r'' -- list of sequences representing r_i^w;
    - ''A'' -- GIM;
    
    Output: ''L'' -- matrix containing the l-vectors:
    '''
    n = len(r)
    l = []
    for i in [0..n-1]:
        l.append(lVector(r[i],n,A))
    
    L = matrix(l)
    
    return L

def betaMutation(betaM, M, w):
    r'''
    Calculates the Beta^w matrix given a mutation sequence, B^w, and C^w
    
    Input:
    - ''betaM'' -- n x n current Beta matrix;
    - ''M'' -- n x 2n current matrix representing labelled seed;
    - ''w'' -- list representing mutation sequence;
    
    Output: ''beta'' -- n x n Beta matrix beta^w:
    '''
    n = M.nrows()
    beta = copy(betaM)
    length = len(w)
    
    if length == 0:
        return beta
    
    for i in w:
        if i < 1 or i > n: # Checks to see if we are mutating at non-vertices
            print("Invalid mutation at vertex: ", i)
            return beta
    
    if length == 1:
        k = w[0]-1
        mut = copy(beta)
        for i in [0..n-1]:
            Bik = M[i,k]
            Ck = M[k,n:]
            BikCk = Bik * Ck
            for j in [0..n-1]:
                Bjk  = M[j,k]
                BjkCk = Bjk * Ck
                
                if i == j:
                    mut[i,j] = beta[i,j]
                elif (BikCk <= 0 and BjkCk <=0) or (BikCk > 0 and BjkCk > 0):
                    mut[i,j] = beta[i,j]
                else:
                    mut[i,j] = beta[i,j] + beta[i,k]*beta[k,j]
                        
        
        beta = copy(mut)
    elif length > 1:
        i = w[length-1]
        v = w[0:length-1]
        MPrime = mMutation(M,v)
        betaPrime = betaMutation(beta,M,v)
        beta = betaMutation(betaPrime, MPrime, [i])
    
    return beta

def betaFromR(A, r):
    r'''
    Calculates the beta matrix from the r_i^w and L^w
    Gives us two ways to calculate the matrix
    
    Input:
    - ''A'' -- GIM;
    - ''r'' -- list of r_i^w's;
    
    Output: ''beta'' -- n x n Beta matrix beta^w:
    '''
    n = len(r)
    L = lMatrix(r,A)
    aM = -(A.transpose())
    lM = L.transpose()
    g = [identity_matrix(n) for i in [0..n-1]]
    for i in [0..n-1]:
        index = floor(len(r[i])/2)
        rI = r[i]
        if index > 0:
            g[i] = (rAction(rI[0:index],A)).inverse()
            
    mut = [[(aM[j,:]*g[j]*lM[:,i])[0][0] for j in [0..n-1]] for i in [0..n-1]]
    beta = matrix(mut)
    
    return beta

def gCancel(rSeq):
    r'''
    Function that gives us g_i^w in the least number of terms but keeps g_i^w well-defined
    
    Input:

    - ''rSeq''  -- sequence of indices; Represents the indices of the product ''r_{i_1}r_{i_2}\\dots r_{i_n}''

    Output: ''red'' -- reduced sequence of indices
    '''
    m = len(rSeq)
    red = rSeq
    if m < 5:
        return red
    
    index = floor(m/2)
    red1 = rCancel(rSeq[0:index])
    red2 = [rSeq[index]]
    red3 = rCancel(rSeq[index+1:m])

    red = red1 + red2 + red3
    
    return red

def rMutation2(i, M, w):
    r'''
    This function takes in an index i, the initial exchange matrix, the initial c-vector matrix, and a mutation sequence w

    Input:
    - ''i'' -- index; Indexing begins at 1 not 0, as we are accustomed to

    - ''M'' -- n x 2n matrix representing labelled seed we are mutating from;
    
    - ''w'' -- mutation sequence;

    Output: ''rW'' -- an index sequence representing ''r_i^w''
    '''
    wPrime = []
    length = len(w)
    rW = [i]
    n = M.nrows()
    
    if length == 0:
        return rW
    
    for j in w:
        if i < 1 or i > n: # Checks to see if we are mutating at non-vertices
            print("Invalid mutation at vertex: ", i)
            return beta
    
    if length == 1:
        k = w[0]-1
        j = i-1
        Bjk = M[j,k]
        Ck = M[k,n:]
        value = Bjk * Ck
        
        if value > 0:
            rW = [k+1, i, k+1]
        else:
            rW = [i]
    else:
        wPrime = w[0:length-1]
        k = w[length-1]-1
        j = i-1
        
        if j == k:
            rW = rMutation2(i,M,wPrime)
        else:
            MPrime = mMutation(M, wPrime)
            #print("M^" + str(wPrime) + ":")
            #print(MPrime)
            Bjk = MPrime[j,k]
            Ck = MPrime[k,n:]
            rW = rMutation2(i,M,wPrime)
            value = Bjk * Ck
            #print("B_" + str(j+1) + str(k+1) +"^"+ str(wPrime) + "C_" + str(k+1) + "^" + str(wPrime) + " = " + str(value))
            if value > 0:
                rW2 = rMutation2(k+1, M, wPrime)
                #print("rW = " + str(rW))
                #print("rW2 = " + str(rW2))
                rW = rW2 + rW + rW2
    #rW = rCancel(rW)
    return rW

'''
Below this point is the positive c vector stuff
'''

def cSigns(M,w):
    r'''
    This function spits out the signs of the c-vectors given an initial framed quiver.
    
    Input:
    - ''M'' -- n x 2n skew-symmetrizable matrix;
    - ''w'' -- mutation sequence
    
    Output:
    - ''c'' -- n-vector signs
    '''
    n = M.nrows()
    C = mMutation(M,w)[:,n:]
    c = [0 for i in [0..n-1]]
    
    for i in [0..n-1]:
        if C[i] > 0:
            c[i] = 1
        elif C[i] < 0:
            c[i] = -1
            
    return c

def spiralPattern(M,w):
    r'''
    This function spits out the spiral pattern of the reflections given an initial framed quiver.
    
    Input:
    - ''M'' -- 3 x 6 skew-symmetrizable matrix;
    - ''w'' -- reduced mutation sequence
    
    Output:
    - ''spirals'' -- n-vector signs
    '''
    n = M.nrows()
    if n != 3:
        print("Not a rank 3 quiver")
        return []
    
    orient = 1
    
    if M[0,1] < 0:
        orient = -1
    
    s = 0
    w = rCancel(w)
    foo = [i for i in [1..n]]
    spirals = [1 for i in [1..n]]
    
    if orient > 0:
        if w[0] == 1:
            s = 2
        elif w[0] ==2:
            s = 3
        elif w[0] == 3:
            s = 1
    if orient < 0:
        if w[0] == 1:
            s = 3
        elif w[0] == 2:
            s = 1
        elif w[0] == 3:
            s = 2
        
    for i in [0..n-1]:
        if rMutation2(i+1,M,w)[0] != s:
            spirals[i] = 0
    
    return spirals

def countEntries(seq, i):
    r'''
    This function counts the number of i's in seq
    
    Input:
    - ''seq'' -- list
    - ''i'' -- integer
    
    Output:
    - ''count'' -- number of times i appears
    '''
    m = len(seq)
    count = 0
    for j in [0..m-1]:
        if seq[j] == i:
            count += 1
    
    return count

def reducedMutations(minLength, maxLength, n):
    r'''
    This function outputs all reduced mutations sequences of certain lengths in a range
    
    Input:
    - ''minLength'' -- the minimum length of our mutation sequences
    - ''maxLength'' -- the maximum length of our mutation sequences
    - ''n'' -- the number of indices
    
    Output:
    - ''mutations'' -- the list containing all of our mutation sequences
    '''
    if minLength > maxLength or minLength < 1:
        minLength = 1
    
    mutations = []
    if maxLength == 1:
        mutations = [[i] for i in [1..n]]
    elif maxLength > 1:
        preMutations = reducedMutations(minLength, maxLength-1,n)
        mutations = preMutations.copy()
        for mut in preMutations:
            if len(mut) == maxLength - 1:
                j = mut[len(mut)-1]
                for i in [1..n]:
                    if i == j:
                        continue
                    else:
                        mutations.append(mut + [i])

        postMutations = mutations.copy()
        if minLength == maxLength:
            for mut in mutations:
                if len(mut) < minLength:
                    postMutations.remove(mut)
                    
        mutations = postMutations
    
    return mutations

def checkLine(rW, orientation, slope, intercept):
    r'''
    This function sees if a line agrees with a certain rW.

    Input:
    - ''rW'' -- the r_i^w sequence with spirals removed
    - ''orientation'' -- the orientation of the universal torus cover + the quiver as a list [o,s,d,v]
    - ''slope'' -- the slope we are supposing that the line has
    - ''intercept'' -- the y-intercept of the supposed line
    
    Output:
    - ''line'' -- true if the line agrees; false otherwise
    '''
    hori = orientation[1]
    diag = orientation[2]
    vert = orientation[3]
    line = true
    if rW[0] != diag:
        print(rW[0])
        print(diag)
        print("The diagram doesn't work because the labelling is off")
        line = false
        return line
    
    if slope == 0: # Remove later
        return line
    
    for i in [1..len(rW)-1]:
        r = rW[:i]
        numH = countEntries(r, hori)+1
        numV = countEntries(r, vert)+1 #These give us the upper right corner of the square we are in

        if rW[i-1] == diag and rW[i] == vert: #If the diagonal comes first, must be in the upper half of the square.
            #Check to see if numH -1 < y < numH is satisfied
            x = numV
            y = slope * x + intercept
            if y >= numH or y <= numH - 1:
                line = false
                print(f"Breaks at {i}; x = {x}; y = {y}; numH = {numH}; numV = {numV}")
                print(r)
                break
        elif rW[i-1] == diag and rW[i] == hori:
            #Check to see if numV -1 < x < numV is satisfied
            y = numH
            x = (y - intercept)/slope
            if x >= numV or x <= numV - 1:
                line = false
                print(f"Breaks at {i}")
                print(r)
                break
        elif rW[i-1] == hori and rW[i] == diag:
            #Check to see if numV - 1 < x < numV is satisfied
            y = numH-1
            x = (y - intercept)/slope
            if x >= numV or x <= numV - 1:
                line = false
                print(f"Breaks at {i}")
                print(r)
                break
        elif rW[i-1] == vert and rW[i] == diag:
            #Check to see if numH - 1 < y < numH is satisfied
            x = numV-1
            y = slope * x + intercept
            if y >= numH or y <= numH - 1:
                line = false
                print(f"Breaks at {i}")
                print(r)
                break
        else:
            line = false
            print("This is not a line")
            return line

    return line

def findLine(rW, orientation, numTrials, fineness):
    r'''
    This function finds a line that agrees with a certain rW.
    If the two lines agree at the midpoint, then we only need the y-intercept to find the line we are trying 

    Input:
    - ''rW'' -- the r_i^w sequence with spirals removed
    - ''orientation'' -- the orientation of the universal torus cover + the quiver as a list [o,s,d,v]
    - ''numTrials'' -- the number of times we try to find a line
    - ''fineness'' -- the denominator of the fractions we try

    Output:
    - ''line'' -- true if the line agrees; false otherwise
    '''
    hori = orientation[1]
    diag = orientation[2]
    vert = orientation[3]
    numHori = countEntries(rW,hori) + 1
    numVert = countEntries(rW,vert) + 1
    m = numHori/numVert
    line = true
    
    if rW[0] != diag:
        print("The diagram doesn't work because the labelling is off")
        line = false
        return line
    
    for j in [0..2*numTrials]:
        y1 = (j-numTrials)/fineness
        slope = m - 2*(y1/numVert)
        line = true
        
        for i in [1..len(rW)-1]:
            r = rW[:i]
            numH = countEntries(r, hori)+1
            numV = countEntries(r, vert)+1 #These give us the upper right corner of the square we are in

            x = 0
            y = 0
            
            if rW[i-1] == diag and rW[i] == vert: #If the diagonal comes first, must be in the upper half of the square.
                #Check to see if numH -1 < y < numH is satisfied
                x = numV
                y = slope * x + y1
                if y >= numH or y <= numH - 1:
                    line = false
                    break
            elif rW[i-1] == diag and rW[i] == hori:
                #Check to see if numV -1 < x < numV is satisfied
                y = numH
                x = (y-y1)/slope
                if x >= numV or x <= numV - 1:
                    line = false
                    break
            elif rW[i-1] == hori and rW[i] == diag:
                #Check to see if numV - 1 < x < numV is satisfied
                y = numH-1
                x = (y-y1)/slope
                if x >= numV or x <= numV - 1:
                    line = false
                    break
            elif rW[i-1] == vert and rW[i] == diag:
                #Check to see if numH - 1 < y < numH is satisfied
                x = numV-1
                y = slope * x + y1
                if y >= numH or y <= numH - 1:
                    line = false
                    break
            else:
                line = false
                print("This is not a line")
                return line
 
        if line:
            print("A line that fits the r_i^w was found: y = " + str(slope) + "x + (" + str(y1) + ")")
            print(checkLine(rW,orientation,slope,y1))
            return line
        
    print("No line was found. Perhaps we should increase the number of trials")
    return line

def findSlopes(rWs, orientation):
    r'''
    This function finds the slopes of all the r_i^w
    
    Input:
    - ''rW'' -- the r_i^w's we are looking at
    - ''orientation'' -- the orientation of the universal torus cover + the quiver as a list [o,s,d,v]
    o stands for positive or negative orientation
    s is the spiral
    d is the diagonal
    v is the vertical
    
    Output:
    - ''slopes'' -- the 3 slopes of r_i^w
    '''
    o = orientation[0]
    s = orientation[1]
    d = orientation[2]
    v = orientation[3]
    denominator = [countEntries(rWs[i],v)+1 for i in [0..2]]
    numerator = [countEntries(rWs[i],s)+1 for i in [0..2]]

    for i in [0..2]:
        if rWs[i][0] == s:
            numerator[i] -= 2
            rWs[i] = rWs[1:len(rWs)-1]
            
    slopes = [[numerator[i], denominator[i]] for i in [0..2]]
    
    return slopes

def slopeDet(slopes):
    r'''
    This function finds the determinant of the 'triangle' formed by the 3 slopes
    
    Input:
    - ''slopes'' -- a list of the form [[a,b], [c,d], [e,f]]
    
    Output:
    - ''det'' -- the determinant of the the 3 x 3 matrix with first column 1's given by the slopes
    '''
    M = matrix([[1] + slopes[i] for i in [0..2]])
    det = M.determinant()
    
    return det

def quiverInfo(M,w):
    r'''
    This function prints out all the collected info that we want
    
    Input:
    - ''M'' -- 3 x 6 skew-symmetrizable matrix;
    - ''w'' -- reduced mutation sequence
    
    Output:
    - test -- true if success or false if not
    '''
    n = M.nrows()
    if n != 3 or len (w) < 2:
        print("Not a rank 3 quiver or mutation sequence too short")
        return false
    
    u = w[:len(w) - 1]
    MPrime = mMutation(M,u)
    B = MPrime[:,:n]
    C = MPrime[:,n:]
    cS = cSigns(M,w)
    sS = spiralPattern(M,w)
    rU = [rCancel(rMutation2(i,M,u)) for i in [1..3]]
    rW = [rCancel(rMutation2(i, M, w)) for i in [1..3]]
    i = w[len(w) -1]
    k = 0
    
    orient = 1
    
    if M[0,1] < 0:
        orient = -1
    
    s = 0 #This determines which of the r_i's count as a spiral/horizontal
    d = 0 #This determines which of the r_i's count as a diagonal
    v = 0 #This determines which of the r_i's count as a vertical
    
    #The following code gives us the correct labelling of the universal cover using the following:
    #The spiral always forms the horizontal components
    #The diagonal is always the first mutation
    #The vertical is the piece left over. As 1,2, & 3 adds up to 6, we can solve for the labelling easily.
    
    if orient > 0:
        s = (1 + w[0]) % 3
    if orient < 0:
        s = (2 + w[0]) % 3
        
    if s == 0:
        s = 3
        
    d = w[0]
    v = 6 - (d + s)
    orientation = [orient,s,d,v]
        
    #Determining the slope:
    #We say that the numerator is the number of r_s plus 1
    #The denominator is the number of r_v plus 1
    
    print("----------")   
    print("Mutation: w = " + str(w))
    print("C-vectors: " + str(cS))
    print("Spirals: " + str(sS))
    print("----------")
    print("Spirals are given by r_" + str(s))
    print("Diagonals are given by r_" + str(d))
    print("Verticals are given by r_" + str(v))
    print("----------")
    
    for j in [0..2]:
        if i - 1 != j and B[j,i-1] * C[i-1,:] > 0:
            k = j+1
            print("r_" + str(k) + "^w = r_" + str(i) + "^u r_" + str(k) + "^u r_" + str(i) + "^u")
            break
    
    slopes = [[1,1],[1,1],[1,1]]
    
    countU = [[countEntries(rU[j], l) for l in [1..3]] for j in [0..2]]
    countW = [[countEntries(rW[j], l) for l in [1..3]] for j in [0..2]]
    print("----------")
    for j in [1..3]:
        print("r_" + str(j) + "^w:")
        print(rW[j-1])
        if j == k:
            print("=")
            print(str(rU[i-1]) + " + " + str(rU[k-1]) + " + " + str(rU[i-1]))
            print("Cancelled portion: ")
            seq = []
            for l in [0..(len(rU[k-1])-3)/2]:
                if rU[k-1][l] == rU[i-1][len(rU[i-1])-1-l]:
                    seq = seq + [rU[k-1][l]]
                else:
                    break
            print(seq)
            if seq != []:
                seqR = list(reversed(seq))
                if seq + [k] + seqR == rU[k-1]:
                    print("The cancelled portion is the entirety of r_" + str(j) + "^u")
                else:
                    print("Conjecture failed.")

        numerator = 0
        denominator = 1
        
        for l in [0..2]:
            print("Number of " + str(l+1) + "'s: " + str(countW[j-1][l]))
        
        print("----------")   
        
        denominator = countW[j-1][v-1]+1
        numerator = countW[j-1][s-1]+1
        r = rW[j-1].copy()
        
        if rW[j-1][0] == s:
            numerator -= 2
            r = r[1:len(r)-1]
        
        
        print("The slope, unsimplified, is given by " + str(numerator) + "/" + str(denominator))
        slopes[j-1] = [numerator,denominator]
        if gcd(numerator, denominator) > 1:
            print("Finding possible line that matches r_" + str(j) + "^w")
            findLine(r, orientation, 10000, 10000)
        elif r != []:
            if not checkLine(r,orientation,numerator/denominator,0):
                findLine(r, orientation, 10000, 10000)
            
        print("----------")  
    print("Slopes: " + str(slopes))
    print("----------")
    slopeM = matrix([[1] + slopes[i] for i in [0..2]])
    print("Determinant of slopes: " + str(slopeM.determinant()))
    print("----------")
    print("Slopes by function: " + str(findSlopes(rW,orientation)))
    print("Determinant by function: " + str(slopeDet(slopes)))
    return true

def matrixR(rSeq):
    '''
    This function calculates the intersection matrix given the reflections
    
    Input:
    - ''rSeq'' -- list of reflections
    
    Output:
    - R -- n x n matrix, where n is the number of reflections
    '''
    n = len(rSeq)
    r = [rCancel(rSeq[i]) for i in [0..n-1]]
    R = matrix([[countEntries(r[i],j) for j in [1..n]] for i in [0..n-1]])
    
    return R


def newQuiverInfo(M,w):
    '''
    This function should represent the recent changes we have made to how we perceive the collected info
    
    Input:
    - ''M'' -- 3 x 6 skew-symmetrizable matrix;
    - ''w'' -- reduced mutation sequence
    
    Output:
    - test -- true if success or false if not
    '''
    n = M.nrows()
    
    u = w[:len(w) - 1]
    MPrime = mMutation(M,w)
    B = MPrime[:,:n]
    C = MPrime[:,n:]
    cS = cSigns(M,w)
    sS = spiralPattern(M,w)
    rU = [rCancel(rMutation2(i,M,u)) for i in [1..n]]
    rW = [rCancel(rMutation2(i, M, w)) for i in [1..n]]
    i = w[len(w) -1]
    k = 0
    
    print("----------")
    print("Matrix M^w:")
    print(MPrime)
    print("Mutation: w = " + str(w))
    print("C-vectors: " + str(cS))
    print("Spirals: " + str(sS))
    print("----------")
    
    for j in [0..n-1]:
        if i - 1 != j and B[j,i-1] * C[i-1,:] > 0:
            k = j+1
            print("r_" + str(k) + "^w = r_" + str(i) + "^u r_" + str(k) + "^u r_" + str(i) + "^u")
            break
    
    R = [[1,1,1],[1,1,1],[1,1,1]]
    
    print("----------")
    for j in [1..n]:
        print("r_" + str(j) + "^w:")
        print(rW[j-1])
        if j == k:
            print("=")
            print(str(rU[i-1]) + " + " + str(rU[k-1]) + " + " + str(rU[i-1]))
            print("Cancelled portion: ")
            seq = []
            for l in [0..(len(rU[k-1])-3)/2]:
                if rU[k-1][l] == rU[i-1][len(rU[i-1])-1-l]:
                    seq = seq + [rU[k-1][l]]
                else:
                    break
            print(seq)

        print("----------")   
   
    R = matrixR(rW)
    print("The interesection matrix R is given by: ")
    print(R)
    return true

def matrixPsi(R, alpha, beta, sigma):
    '''
    This function should produce the Psi Matrix corresponding to an R
    
    Input:
    - ''R'' -- n x n intersection matrix;
    - ''alpha'' -- list of alpha's;
    - ''beta'' -- list of beta's
    - ''sigma'' -- 0 or 1
    
    Output:
    - M -- the matrix Psi
    '''
    n = R.nrows()

    maxIndex = [0 for i in [0..n-1]]
    maxValue = [0 for i in [0..n-1]]
    M = copy(R)
                
    for i in [0..n-1]:
        for j in [0..n-1]:
            M[i,j] = M[i,j]+1
    M[alpha[n-1],beta[0]] = M[alpha[0],beta[0]] + M[alpha[1],beta[0]]
    M[alpha[n-1],beta[1]] = M[alpha[0],beta[1]] + M[alpha[1],beta[1]] - 2*(1-sigma)
    for i in [0..n-1]:
        M[alpha[i],beta[n-1]] = 1

    return M

def intersectionMutation(rList, w):
    r'''
    This function produces the intersection tuple formed by mutating along the mutation sequence w
    
    Input:
    - ''rList'' -- 3 x 3 intersection matrix, 3x1 list of signs, positive or negative sigma
    - ''w'' -- mutation sequence
    
    Output:
    - ''rMut'' -- mutated r-tuple
    '''
    length = len(w)
    rMut = copy(rList)
    
    if length <= 0:
        return rMut
    elif length == 1:
        R = copy(rMut[0])
        RPrime = copy(rMut[0])
        chis = copy(rMut[1])
        chisPrime = copy(rMut[1])
        sigma = rMut[2]
        k = w[0]-1
        for i in [0..2]:
            if i == k:
                chisPrime[i] = -1*chis[i]
                continue
            for j in [0..2]:
                if -1*sigma*chis[k]*R[k,j] > 0 and ((i == 0 and k == 2) or (i == 1 and k == 0) or (i == 2 and k == 1)):
                    value = chis[i]*R[i,j] + 2*chis[k]*R[k,j]
                    RPrime[i,j] = abs(value)
                    if RPrime[i,j] > 0:
                        chisPrime[i] = RPrime[i,j]/value
                elif sigma*chis[k]*R[k,j] > 0 and ((i == 0 and k == 1) or (i == 1 and k == 2) or (i == 2 and k == 0)):
                    value = chis[i]*R[i,j] + 2*chis[k]*R[k,j]
                    RPrime[i,j] = abs(value)
                    if RPrime[i,j] > 0:
                        chisPrime[i] = RPrime[i,j]/value
                else:
                    continue
        rMut = [RPrime,chisPrime, -1*sigma]
    else:
        rMut = intersectionMutation(intersectionMutation(rList,w[:length-1]), [w[length-1]])
        
    return rMut