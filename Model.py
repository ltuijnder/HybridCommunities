# this scripts holds the necessary function to create the matrix elements you want.
import numpy as np

def IsStable(JacobianMatrix):
    EigVals=np.linalg.eigvals(JacobianMatrix)
    EigValsReal=EigVals.real # Only work with real elements
    if(EigValsReal[EigValsReal>0].size==0): # No elements are greater than 0
        return 1
    else:
        return 0

def Create_a(n,Cup):
    Clo=np.transpose(Cup)# Lower triangle connection (i>j)
    C=Cup+Clo# Symmetric connection
    #Create some array slicing, and rightaway also determine it's size
    Upper_r=Cup==-1;sUr=Upper_r[Upper_r==True].size
    Lower_s=Clo==-2;sLs=Lower_s[Lower_s==True].size
    Upper_t=Cup==-3;sUt=Upper_t[Upper_t==True].size
    Lower_t=Clo==-3;sLt=Lower_t[Lower_t==True].size
    m=C==1;sm=m[m==True].size
    
    # Create the random Coëfficients. 
    e=np.random.rand(sm)
    gt=np.random.rand(sLt)
    gs=np.random.rand(sLs)
    
    A=np.random.rand(n,n)
    A[C==0]=0 # Keep only the non negatif ones.
    ASum=np.sum(A,1)
    
    row_ind=np.indices((n,n))[0]# Row indices. We need this since we will use this to acces ASum. 
    # Also this is in nxn matrix which needs to be, since we need to array slices is, which are also nxn matrixes.
    # the result will be an array with the correspoding row indices of the elements. Then apply this to ASum. 
    
    a=np.zeros((n,n))
    #First antagonistic
    # t
    a[Upper_t]=A[Upper_t]/ASum[row_ind[Upper_t]]
    a+=np.transpose(a)# Symetrise it. Because of this step, we need t first, since now it is still empty.
    a[Lower_t]*=gt # This should acces the term we just transposed there.
    a[Upper_t]*=-1 #Still make the upper negatif
    # r
    a[Upper_r]=-1*A[Upper_r]/ASum[row_ind[Upper_r]]
    # s
    a[Lower_s]=gs*A[Lower_s]/ASum[row_ind[Lower_s]]
    
    # Second Mutalistic reactions
    a[m]=e*A[m]/ASum[row_ind[m]]
    
    return a

def Generate(n,P,fm,fr,fs):# see other notebook for comments.
    p_arr=np.array([P*(1-fm)*fr,P*(1-fm)*fs,P*(1-fm)*(1-fr-fs),P*fm,1-P])# Chance distribution fo the connection type
    v=np.array([-1,-2,-3,1,0])# [r,s,t,m,nothing] This is the code that we will use in the C matrix.
    Pair_type=np.random.choice(v,int(n*(n-1)/2),p=p_arr)
    Cup=np.zeros((n,n))
    Cup[np.triu_indices(n,1)]=Pair_type# upper triangle aka (i<j) (negatif for antagonistic)
    
    a=Create_a(n,Cup)
    
    X_fixed=np.random.rand(n)
    s=np.random.rand(n)
    
    J=np.dot(np.diag(X_fixed),a)+np.diag(-1*s*X_fixed)# first random is saturation, second is x_fixed
    StableBool=IsStable(J)
    
    # complete a fully:
    a+=np.diag(-1*s)
    
    #compute r
    r=-1*np.dot(a,X_fixed)
    
    return (r,a,X_fixed,StableBool)
        