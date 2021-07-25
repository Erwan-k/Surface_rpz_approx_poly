
from math import sin,cos,pi,atan2,exp
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

def rotation_d_angle_teta_autour_de_U(Vect,U,teta):
    U = normer(U)
    V = rechercher_V_orthogonal_à_(U)
    V = normer(V)
    W = produit_vectoriel(U,V)
    P = [[U[0],V[0],W[0]],[U[1],V[1],W[1]],[U[2],V[2],W[2]]]
    Pinv = [U,V,W]
    M = [[1,0,0],[0,cos(teta),-sin(teta)],[0,sin(teta),cos(teta)]]
    R = produit_matriciel(P,produit_matriciel(M,Pinv))
    return produit_matrice_vecteur(R,Vect)

def normer(U): #s = (sum([i**2 for i in U]))**(1/2)
    s = 0
    for i in U:
        s+=i**2
    s = s**(1/2)
    return [U[0]/s,U[1]/s,U[2]/s]



def rechercher_V_orthogonal_à_(U):
    if U == [1,0,0] or U == [-1,0,0] or U == [0,1,0] or U == [0,-1,0]:
        V = [0,0,1]
    elif U == [0,0,1] or U == [0,0,-1]:
        V = [1,0,0]
    else:
        if U[0] == 0:
            V = [0,-U[2],U[1]]
        else:
            a = -U[2]
            c = U[0]
            V = [a,0,c]
    return V

def produit_vectoriel(U,V):
    [xu,yu,zu] = U
    [xv,yv,zv] = V
    X = yu*zv-yv*zu
    Y = zu*xv-zv*xu
    Z = xu*yv-xv*yu
    return [X,Y,Z]

def produit_matriciel(A,B):
    [[A11,A12,A13],[A21,A22,A23],[A31,A32,A33]] = A
    [[B11,B12,B13],[B21,B22,B23],[B31,B32,B33]] = B
    C11 = A11*B11+A12*B21+A13*B31
    C12 = A11*B12+A12*B22+A13*B32
    C13 = A11*B13+A12*B23+A13*B33
    
    C21 = A21*B11+A22*B21+A23*B31
    C22 = A21*B12+A22*B22+A23*B32
    C23 = A21*B13+A22*B23+A23*B33
    
    C31 = A31*B11+A32*B21+A33*B31
    C32 = A31*B12+A32*B22+A33*B32
    C33 = A31*B13+A32*B23+A33*B33

    return [[C11,C12,C13],[C21,C22,C23],[C31,C32,C33]]

def produit_matrice_vecteur(A,B):
    [[A11,A12,A13],[A21,A22,A23],[A31,A32,A33]] = A
    [x,y,z] = B
    X = A11*x+A12*y+A13*z
    Y = A21*x+A22*y+A23*z
    Z = A31*x+A32*y+A33*z
    return [X,Y,Z]


#####

def projection_sur_P_parallellement_a_D(n,u,pt):
    [a,b,c],[xu,yu,zu],[x,y,z] = n,u,pt
    
    s = a*xu+b*yu+c*zu
    xp = x-(x*a*xu/s+y*b*xu/s+z*c*xu/s)
    yp = y-(x*a*yu/s+y*b*yu/s+z*c*yu/s)
    zp = z-(x*a*zu/s+y*b*zu/s+z*c*zu/s)

    return [xp,yp,zp]

def projection_ortho_sur_P_dirige_par_u_(u,pt):
    return projection_sur_P_parallellement_a_D(u,u,pt)

def projection_ortho_sur_oyz(pt):
    return projection_ortho_sur_P_dirige_par_u_([1,0,0],pt)


#####
    

def afficher_segments_vu_selon_la_direction_u(segment,u,couleur,epaisseur):
    seg = deepcopy(segment)
    
    #Je cherche un vecteur ortho_n orthogonal a u et la projection de u sur Oxy.
    ortho_n = produit_vectoriel(n,[n[0],n[1],0])
    
    #Je récupère les coordonnées sphériques de n
    alpha,beta = atan2(n[1],n[0]),atan2(n[2],(n[0]**2+n[1]**2)**(1/2))

    #A chaque point : 
    for i in range(len(seg)):
        for j in range(len(seg[i])):
            
            #Si beta est non nul
            if beta:
                #Je fais une rotation d'angle -beta autour de ortho_n 
                seg[i][j] = rotation_d_angle_teta_autour_de_U(seg[i][j],ortho_n,-beta)
            
            #Si alpha est non nul
            if alpha:
                #Je fais une rotation d'angle -alpha autour de Uz
                seg[i][j] = rotation_d_angle_teta_autour_de_U(seg[i][j],[0,0,1],-alpha)
            
            #Je projette le point sur Oyz, parallèlement à 
            seg[i][j] = projection_ortho_sur_oyz(seg[i][j])
        
            #Je ne conserve que les coordonnées y et z de la projection
            seg[i][j] = [seg[i][j][1],seg[i][j][2]]

    #J'affiche les segments
    for i in seg:
        plt.plot([i[0][0],i[1][0]],[i[0][1],i[1][1]],couleur,linewidth=epaisseur)

    
        
    


def f(x,y):
    return (x**2-y**2)*exp(-(x**2+y**2))

def g(x,y):
    #return -exp(-1)+x**2+(y-1)**2 #approximation en (0;1)
    return -exp(-1)+x**2+(y+1)**2 #approximation en (0;1)
    #return x**2-y**2 #approximation en (0;0)
    #return exp(-1)-(x+1)**2-y**2 #approximation en (-1;0)
    #return exp(-1)-(x-1)**2-y**2 #approximation en (-1;0)

seg3 = []
etendue_x,etendue_y = 2.5,2.5
X,Y = np.linspace(-etendue_x,etendue_x,31), np.linspace(-etendue_y,etendue_y,31)
for i in range(len(X)-1):
    for j in range(len(Y)):
        seg3 += [[  [X[i],Y[j],f(X[i],Y[j])] , [X[i+1],Y[j],f(X[i+1],Y[j])]  ]]
for i in range(len(X)):
    for j in range(len(Y)-1):
        seg3 += [[  [X[i],Y[j],f(X[i],Y[j])] , [X[i],Y[j+1],f(X[i],Y[j+1])]  ]]
        
seg4 = []
a,b = 0,-1
etendue_x,etendue_y = 1,1
X,Y = np.linspace(a-etendue_x,a+etendue_x,16), np.linspace(b-etendue_y,b+etendue_y,16)
for i in range(len(X)-1):
    for j in range(len(Y)):
        seg4 += [[  [X[i],Y[j],g(X[i],Y[j])] , [X[i+1],Y[j],g(X[i+1],Y[j])]  ]]
for i in range(len(X)):
    for j in range(len(Y)-1):
        seg4 += [[  [X[i],Y[j],g(X[i],Y[j])] , [X[i],Y[j+1],g(X[i],Y[j+1])]  ]]


#Je veux tracer ces segments
seg1 = [[[-10**6,0,0],[10**6,0,0]],
        [[0,-10**6,0],[0,10**6,0]],
        [[0,0,-10**6],[0,0,10**6]]]
       
alphas = np.linspace(0,2*pi,51)
zs = [sin(12*pi*i/len(alphas)) for i in range(len(alphas))]
for i in range(len(alphas)):

    n = [3*cos(alphas[i]),3*sin(alphas[i]),0.3]
       
    plt.axis([-3,3,-0.7,0.7])

    afficher_segments_vu_selon_la_direction_u(seg1,n,"k",0.5)
    afficher_segments_vu_selon_la_direction_u(seg4,n,"r",1)
    afficher_segments_vu_selon_la_direction_u(seg3,n,"b",1)
    
    num_image = str(i)
    num_image = "0"*(3-len(num_image))+num_image
    plt.savefig("image"+num_image+".png")
    print(num_image)
    
    plt.show()    





























