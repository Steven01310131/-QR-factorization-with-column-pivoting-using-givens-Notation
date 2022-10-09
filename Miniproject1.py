
###############################################################################################
# Miniproject 1 Scientific computing bridging cource                                          #
# Author: Stefanos Tsampanakis, stefanos.tsabanakis@gmail.com                                 #
# Python function for QR factorization with column pivoting using givens Notation:            #
# input matrix A(m x n) with m>=n and output Q,R,P, and r, where r is the numerical rank of A #
###############################################################################################
from xml.dom.minidom import Notation
import numpy as np
from scipy import linalg
import scipy
import rogues 
A = np.array([[1,-1,2,0],[1,2,-1,3],[1,1,0,2],[1,-1,2,0],[1,3,-1,4]])
def QRcPivot(A):
    (num_rows,num_cols) = np.shape(A)
    P=[]
    G=np.identity(num_rows)
    for n in range(num_cols):
        k=0
        I=np.identity(num_cols)
        II=I[n:,n:]#submatrix of identity
        sum3=0
        B=A[n:,n:]#submatrix of A
        (sub_num_rows, sub_num_cols) = np.shape(B)
        for j in range(sub_num_cols) :          #iterative method for finding the collumn with max norm2 and remembering 
            sum1=0                              #in which collumn the max norm is in the submatrix 
            sum2=0                              #               
            for i in range(sub_num_rows):       #   
                sum1=sum1+abs(B[i,j])**2        #   
                sum2=np.sqrt(sum1)              #
                if sum2>sum3 :                  #  
                    sum3=sum2                   #
                    k=j                         #
        II[:, [0, k]] = II[:, [k, 0]]
        I[n:,n:]=II # permutation matrix 
        P=P+[I]
        A=np.dot(A,I) 
        B[:, [0, k]] = B[:, [k, 0]]            
        G2=np.identity(sub_num_rows)
        G3=np.identity(num_rows)
        for j in range(sub_num_rows-1,0,-1):
            c,s=givensparameters(B[0,0], B[j,0])
            G1=np.identity(sub_num_rows)
            G1[0,0]=c
            G1[j,j]=c
            G1[j,0]=s
            G1[0,j]=-s
            G2=np.dot(G2,G1)
            B=np.dot(G1,B)
        G3[n:,n:]=G2

        G=np.dot(G,G3)
        A[n:,n:]=B
    rank=np.linalg.matrix_rank(A)
    return A ,  P,   rank,   G
def givensparameters(a, b):
    r = np.sqrt(a**2+b**2)
    c = a/r
    s = -b/r

    return (c, s)
A = np.array([[1,-1,2,0],[1,2,-1,3],[1,1,0,2],[1,-1,2,0],[1,3,-1,4]])
(A ,  P,   rank,   Q)=QRcPivot(A)
print(A)
#using scipy library 
def QRcPivotlibrary(A):
    (Q,R,P)=scipy.linalg.qr(A,pivoting=True)
    rank=np.linalg.matrix_rank(A)
    return Q,R,P,rank

######################################################################
# Function to solve the least squares problem with pivot             #
######################################################################




def QRcPivotlibrary(A):
    (Q,R,P)=scipy.linalg.qr(A,pivoting=True)
    rank=np.linalg.matrix_rank(A)
    return Q,R,P,rank

def LeastSquaresQRcPivot(A,b):
    (Q,R,P,rank)=QRcPivotlibrary(A)
    Q1=Q[:,:rank]
    R1=R[:rank,:rank]
    L=np.matmul(Q1,b[:rank,:])
    y=np.zeros(rank)
    x=np.zeros(len(b))
    for i in range(rank-1,-1,-1):#backwards
        y[i]=(L[i]-np.sum(R1[i,i+1:rank]*y[i+1:rank]))/R[i,i]
    x[:rank]=y
    for i in P:
        x[[0,i]] = x[[i, 0]]
    return x

A = np.array([[1,-1,2,0],[1,2,-1,3],[1,1,0,2],[1,-1,2,0],[1,3,-1,4]])   



A = np.array([[1,-1,2,0],[1,2,-1,3],[1,1,0,2],[1,-1,2,0],[1,3,-1,4]])
b=np.array([[1],[-1],[0],[1]])
(Q,R,P,rank)=QRcPivotlibrary(A)
print(R)
# Z=np.identity(5)
# # print(Z)
# J=[]
# for i in P:
#     Z=np.identity(4)
    
#     Z[:, [i, 0]] = Z[:, [0, i]]
#     J=J+[Z]
# print(A)
# for i in J:
#     A=np.matmul(A,i)
# print(A-np.matmul(Q,R))
# LeastSquaresQRcPivot(A,b)
n = 100
B = rogues.neumann(n)
A = B[0]

