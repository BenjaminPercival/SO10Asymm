#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 19 12:55:29 2022

@author: wmnc67
"""


from z3 import * 

#import itertools 

NAHE=[[1,1,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0],\
      [1,0,1,0,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1,0],\
      [1,0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,1,1,1,1,1,0,0,1]]

X = [ [ Int("x_%s_%s" % (i+1, j+1)) for j in range(36) ] \
      for i in range(3) ]

RealBCs  = [ Or( X[i][j]==0, X[i][j]==1) \
             for i in range(3) for j in range(36) ]

#supercurrent
#-------------------------------------------------------------------------
Scurrent1 = [(X[i][0]+X[i][1]+X[i][4]+X[i][10])%2==0 \
            for i in range(3) ]
Scurrent2 = [(X[i][0]+X[i][1]+X[i][5]+X[i][11])%2==0 \
            for i in range(3) ]
Scurrent3 = [(X[i][0]+X[i][2]+X[i][6]+X[i][12])%2==0 \
            for i in range(3) ]
Scurrent4 = [(X[i][0]+X[i][2]+X[i][7]+X[i][13])%2==0 \
            for i in range(3) ]
Scurrent5 = [(X[i][0]+X[i][3]+X[i][8]+X[i][14])%2==0 \
            for i in range(3) ]
Scurrent6 = [(X[i][0]+X[i][3]+X[i][9]+X[i][15])%2==0 \
            for i in range(3) ]

Scurrent=Scurrent1+Scurrent2+Scurrent3+Scurrent4+Scurrent5+Scurrent6
#psi,chis 
#-------------------------------------------------------------------------
ChiPsi = [Or(And(X[i][0]==1,X[i][1]+X[i][2]+X[i][3]==1),\
               And(X[i][0]==0,X[i][1]+X[i][2]+X[i][3]==2),\
               X[i][0]+X[i][1]+X[i][2]+X[i][3]==0) \
                   for i in range(3) ]
    

#-------------------------------------------------------------------------
#SO(10)=> psi^12345 BCs same and conistent U(1) charges
Psi12345 = [Or(And(X[i][28+m]==1,X[i][33]+X[i][34]+X[i][35]==1),\
               And(X[i][28+m]==0,X[i][33]+X[i][34]+X[i][35]==2)) \
                   for i in range(2) for m in range(5)]
Phisb61=[sum(X[2][28:32])==4]
Phisb62=[sum(X[2][32:])==0]
b6z1=Phisb61+Phisb62
#-------------------------------------------------------------------------
#dot prods MI with NAHE
NAHEbj=[]
for m in range(3):
    for n in range(2):
        ls=[NAHE[m][j]*X[n][j] for j in range(36)]
        DPmn=(2*sum(ls[:4])+sum(ls[4:16])-sum(ls[16:28])-2*sum(ls[28:]))/2
        NAHEbj.append(DPmn%2==0)
#b6 DPs
for m in range(3):
    lsb6=[NAHE[m][j]*X[2][j] for j in range(36)]
    DPmnb6=(2*sum(lsb6[:4])+sum(lsb6[4:16])-sum(lsb6[16:28]))/2
    NAHEbj.append(DPmnb6%2==0)       

#-------------------------------------------------------------------------
#dot prods MI with each other
b4b5=[]
Prod45=[X[0][j]*X[1][j] for j in range(36)]
DP45=(2*sum(Prod45[:4])+sum(Prod45[4:16])-sum(Prod45[16:28])-2*sum(Prod45[28:]))/2
b4b5.append(DP45%2==0)
#set_option(max_args=10000000, max_lines=1000000, max_depth=10000000, max_visited=1000000)
#print(b4b5)
b4b6=[]
Prod46=[X[0][j]*X[2][j] for j in range(28)]
DP46=(2*sum(Prod46[:4])+sum(Prod46[4:16])-sum(Prod46[16:28]))/2
b4b6.append(DP46%2==0)

b5b6=[]
Prod56=[X[1][j]*X[2][j] for j in range(28)]
DP56=(2*sum(Prod56[:4])+sum(Prod56[4:16])-sum(Prod56[16:28]))/2
b5b6.append(DP56%2==0)

#-------------------------------------------------------------------------
#dot prods MI with themselves
b4b4=[]
Prod44=[X[0][j]*X[0][j] for j in range(36)]
DP44=(2*sum(Prod44[:4])+sum(Prod44[4:16])-sum(Prod44[16:28])-2*sum(Prod44[28:]))/2
b4b4.append(DP44%4==0)

b5b5=[]
Prod55=[X[1][j]*X[1][j] for j in range(36)]
DP55=(2*sum(Prod55[:4])+sum(Prod55[4:16])-sum(Prod55[16:28])-2*sum(Prod55[28:]))/2
b5b5.append(DP55%4==0)

b6b6=[]
Prod66=[X[2][j]*X[2][j] for j in range(36)]
DP66=(2*sum(Prod66[:4])+sum(Prod66[4:16])-sum(Prod66[16:28])-2*sum(Prod66[28:]))/2
b6b6.append(DP66%4==0)

dotProdsMI=NAHEbj+b4b5+b4b6+b5b6+b4b4+b5b5+b6b6
#-------------------------------------------------------------------------
#no real pairings
#b1 group
bkb1L=[[X[i][6+k] for k in range(4)]\
      for i in range(3)]
bkb1R=[[X[i][18+k] for k in range(4)]\
      for i in range(3)]
  
Pairs1L=[sum(ls) for ls in bkb1L] #ensure we get even number on left and right...
Pairs1R=[sum(ls) for ls in bkb1R]

Pair1Cs1L=[x%2==0 for x in Pairs1L]
Pair1Cs1R=[x%2==0 for x in Pairs1R]
Pair1Cs2=[Not(And(Pairs1L[i]==4,Pairs1R[i]==0)) for i in range(3)] 
Pair1Cs3=[Not(And(Pairs1L[i]==0,Pairs1R[i]==4)) for i in range(3)]

bkb1=[[X[i][6+k] for k in range(4)]+[X[i][18+k] for k in range(4)]\
      for i in range(3)]  
bkb1Prods=[[bkb1[m][j]*bkb1[n][j] for j in range(4) if m<n] \
           for m in range(3) for n in range(3) ]

bkb1P2 = [x for x in bkb1Prods if x != []]

Pairs2=[sum(ls)%2==0 for ls in bkb1P2]

Group1=Pair1Cs1L+Pair1Cs1R+Pair1Cs2+Pair1Cs3+Pairs2
#b2 group
bkb2L=[[X[i][4+k] for k in range(2)]+[X[i][14+k] for k in range(2)]\
      for i in range(3)]
bkb2R=[[X[i][16+k] for k in range(2)]+[X[i][26+k] for k in range(2)]\
      for i in range(3)]
Pairs2L=[sum(ls) for ls in bkb2L] #ensure we get even number on left and right...
Pairs2R=[sum(ls) for ls in bkb2R]

Pair2Cs1L=[x%2==0 for x in Pairs2L]
Pair2Cs1R=[x%2==0 for x in Pairs2R]
Pair2Cs2=[Not(And(Pairs2L[i]==4,Pairs2R[i]==0)) for i in range(3)]
Pair2Cs3=[Not(And(Pairs2L[i]==0,Pairs2R[i]==4)) for i in range(3)]

bkb2=[[X[i][4+k] for k in range(2)]+[X[i][14+k] for k in range(2)]+\
      [X[i][16+k] for k in range(2)]+[X[i][26+k] for k in range(2)]\
      for i in range(3)]

bkb2Prods=[[bkb2[m][j]*bkb2[n][j] for j in range(8) if m<n] \
           for m in range(3) for n in range(3) ]
bkb2P2 = [x for x in bkb2Prods if x != []]

Pairs4=[sum(lst)%2==0 for lst in bkb2P2]

Group2=Pair2Cs1L+Pair2Cs1R+Pair2Cs2+Pair2Cs3+Pairs4
#b3 group
bkb3L=[[X[i][10+k] for k in range(4)]\
      for i in range(3)]
bkb3R=[[X[i][22+k] for k in range(4)]\
      for i in range(3)]   
Pairs3L=[sum(ls) for ls in bkb3L]
Pairs3R=[sum(ls) for ls in bkb3R]

Pair3Cs1L=[x%2==0 for x in Pairs3L]
Pair3Cs1R=[x%2==0 for x in Pairs3R]
Pair3Cs2=[Not(And(Pairs3L[i]==4,Pairs3R[i]==0)) for i in range(3)]
Pair3Cs3=[Not(And(Pairs3L[i]==0,Pairs3R[i]==4)) for i in range(3)]

bkb3=[[X[i][10+k] for k in range(4)]+[X[i][22+k] for k in range(4)]\
      for i in range(3)]

bkb3Prods=[[bkb3[m][j]*bkb3[n][j] for j in range(8) if m<n] \
           for m in range(3) for n in range(3) ]
bkb3P2 = [x for x in bkb3Prods if x != []]

Pairs6=[sum(ls)%2==0 for ls in bkb3P2]

Group3=Pair3Cs1L+Pair3Cs1R+Pair3Cs2+Pair3Cs3+Pairs6

PairingConstraints=Group1+Group2+Group3
#-------------------------------------------------------------------------
#Not the same internal fermions in 2
Com_b4b5=[X[0][i]==X[1][i] for i in range(4,28)]
Different1=[Not(And(Com_b4b5[0],Com_b4b5[1],Com_b4b5[2],Com_b4b5[3],Com_b4b5[4],\
                    Com_b4b5[5],Com_b4b5[6],Com_b4b5[7],Com_b4b5[8],Com_b4b5[9],\
                    Com_b4b5[10],Com_b4b5[11],Com_b4b5[12],Com_b4b5[13],Com_b4b5[14],\
                    Com_b4b5[15],Com_b4b5[16],Com_b4b5[17],Com_b4b5[18],Com_b4b5[19],\
                    Com_b4b5[20],Com_b4b5[21],Com_b4b5[22],Com_b4b5[23]))]
Com_b4b6=[X[0][i]==X[2][i] for i in range(4,28)]
Different2=[Not(And(Com_b4b6[0],Com_b4b6[1],Com_b4b6[2],Com_b4b6[3],Com_b4b6[4],\
                    Com_b4b6[5],Com_b4b6[6],Com_b4b6[7],Com_b4b6[8],Com_b4b6[9],\
                    Com_b4b6[10],Com_b4b6[11],Com_b4b6[12],Com_b4b6[13],Com_b4b6[14],\
                    Com_b4b6[15],Com_b4b6[16],Com_b4b6[17],Com_b4b6[18],Com_b4b6[19],\
                    Com_b4b6[20],Com_b4b6[21],Com_b4b6[22],Com_b4b6[23]))]
Com_b5b6=[X[1][i]==X[2][i] for i in range(4,28)]
Different3=[Not(And(Com_b5b6[0],Com_b5b6[1],Com_b5b6[2],Com_b5b6[3],Com_b5b6[4],\
                    Com_b5b6[5],Com_b5b6[6],Com_b5b6[7],Com_b5b6[8],Com_b5b6[9],\
                    Com_b5b6[10],Com_b5b6[11],Com_b5b6[12],Com_b5b6[13],Com_b5b6[14],\
                    Com_b5b6[15],Com_b5b6[16],Com_b5b6[17],Com_b5b6[18],Com_b5b6[19],\
                    Com_b5b6[20],Com_b5b6[21],Com_b5b6[22],Com_b5b6[23]))]
Different=Different1+Different2+Different3
#-------------------------------------------------------------------------
#preserve at least one of first torus moduli, WLOG take h12 
h12=[X[i][4]+X[i][10]+X[i][17]+X[i][23]%2==0 for i in range(3)]


#-------------------------------------------------------------------------
#project the 2nd and 3rd tori moduli
h33=[X[i][6]+X[i][12]+X[i][18]+X[i][24]%2==1 for i in range(3)]
h34=[X[i][6]+X[i][12]+X[i][19]+X[i][25]%2==1 for i in range(3)]
h43=[X[i][7]+X[i][13]+X[i][18]+X[i][24]%2==1 for i in range(3)]
h44=[X[i][7]+X[i][13]+X[i][19]+X[i][25]%2==1 for i in range(3)]
h55=[X[i][8]+X[i][14]+X[i][20]+X[i][26]%2==1 for i in range(3)]
h56=[X[i][8]+X[i][14]+X[i][21]+X[i][27]%2==1 for i in range(3)]
h65=[X[i][9]+X[i][15]+X[i][20]+X[i][26]%2==1 for i in range(3)]
h66=[X[i][9]+X[i][15]+X[i][21]+X[i][27]%2==1 for i in range(3)]

ModuliConstraints=h12+h33+h34+h43+h44+h55+h56+h65+h66
#-------------------------------------------------------------------------
constraints=RealBCs+Scurrent+ChiPsi+Psi12345+b6z1+dotProdsMI\
    +PairingConstraints+Different\
        +ModuliConstraints
#-------------------------------------------------------------------------

s=Solver()
countr=0
s.add(constraints)
if s.check() == unsat:
    print("failed to solve")
else:
    while s.check() == sat:
        m = s.model()
        r = [ [ m.evaluate(X[i][j]) for j in range(36) ] 
              for i in range(3) ]
        #print(r)
        countr+=1
        f = open('BasesB.txt','a') 
    
        old_stdout = sys.stdout  #  store the default system handler to be able to restore it 
    
        sys.stdout = f 
        print(countr)
        for ls in r:
            
            print(ls)
        f.close()
        sys.stdout=old_stdout
        
        s.add(Not(And([v() == m[v] for v in m]))) 

print(countr)