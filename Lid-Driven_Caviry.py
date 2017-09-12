
# Import library
import matplotlib.pyplot as plt
import numpy as np
import enum
import shutil  
import os
# Types of nodes as follows:
class NodeType(enum.IntEnum):

    INTERIOR = 0
    WALL = 1
    LID = 2

'''make geometry'''
def makeGeometry(): 

    # build rectangular domain
    pos_x, pos_y = np.meshgrid(np.linspace(0, Length, MB), 
                               np.linspace(0,  Height, NB))
    
    pos_x = np.transpose(pos_x)
    pos_y = np.transpose(pos_y)

    # give type of all nodes
    nodes_type = np.zeros_like(pos_x)

    for i in range( 0, MB ):
        for j in range( 0, NB ):
            
            if ( pos_x[i][j] != 0 and pos_x[i][j] != Length and
                 pos_y[i][j] != 0 and pos_y[i][j] != Height):
                 
                 nodes_type[ i, j ] = NodeType.INTERIOR

            elif ( pos_y[i][j] == Height):
                 
                 nodes_type[ i, j ] = NodeType.LID 

            else:
                 nodes_type[i][j] = NodeType.WALL

    return ( pos_x, pos_y, nodes_type )


'''define Initail conditions of Psi'''
def IntPsi():
    #assume Psi = 0 at time = 0
    Psi = np.zeros_like(pos_x)
    #assume Vorticity, u, v = 0 at time = 0
    Vor = np.zeros_like(pos_x)
    u   = np.zeros_like(pos_x)
    v   = np.zeros_like(pos_x)

    #assume Vorticity -2*u0/h on lid at time = 0;
    # Vor[:,NB-1] = -2*u0/dx

    return ( Psi, Vor ,u ,v ) 

'''compute Psi'''
def computePsi():

    # interior part
    for i in range(101): 
        Psi[1:MB-1, 1:NB-1] = 0.25*( Vor[1:MB-1, 1:NB-1]*dx*dx + 
                                        (Psi[1:MB-1,0:NB-2] + Psi[1:MB-1,2:NB]) +
                                        (Psi[0:MB-2,1:NB-1] + Psi[2:MB,1:NB-1])   )

'''compute Vel from Psi'''        
def computeVel():
    # v = -dPsi/dx
    for i in range (1,MB):
        for j in range (NB):
            # skip over walls, otherwise differencing on neighbors will be off
            if (i == 0):
                 v[i,j] = - (Psi[i+1,j]-Psi[i,j]) / dx
            elif (i == MB-1):
                 v[i,j] = - (Psi[i,j]-Psi[i-1,j]) / dx
            else:
                 v[i,j] = - (Psi[i+1,j]-Psi[i-1,j]) / (2*dx)
            
    #u = dpsi/dy
    for i in range (1,MB):
        for j in range (NB):
            if (j == 0):
                 u[i,j] = (Psi[i,j+1] - Psi[i,j]) / dy
            elif ( j == NB-1):   
                 u[i,j] = (Psi[i,j] - Psi[i,j-1]) / dy
            else:
                 u[i,j] = (Psi[i,j+1] - Psi[i,j-1]) / (2*dy)


'''compute RHS for Vorticity equation''' 
def RHS(Vor):
    dx2 = dx*dx
    dy2 = dy*dy    
    # make copy so we use consistent data
    r = np.zeros_like(Vor)

    for i in range(1,MB-1):
        for j in range(1,NB-1):
                
            #viscous term, d^2w/dz^2+d^2w/dr^2+(1/r)dw/dr
            A = nu*(
                    (Vor[i-1, j] - 2*Vor[i, j] + Vor[i+1, j])/dx2 + 
                    (Vor[i, j-1] - 2*Vor[i, j] + Vor[i, j+1])/dy2 )
            
            #convective term u*dw/dx    
            B = - u[i, j]*(Vor[i+1, j] - Vor[i-1, j])/(2*dx)

            #convective term v*dw/dy
            C = v[i, j]*(Vor[i, j+1] - Vor[i, j-1])/(2*dy)
                      
            r[i, j] = A + B + C
    return r

'''compute Vorticity on the boundary''' 
def VorBou():
    dx2 = dx*dx
    dy2 = dy*dy
    # wall boundaries
    for i in range(MB):
        for j in range(NB):

            #left wall
            if ( i == 0 and nodes_type[i,j] == NodeType.WALL):
                Vor[i,j] = 2 * ( Psi[i,j] - Psi[i+1,j] )/ dx2

            #right wall
            if ( i == MB-1 and nodes_type[i,j] == NodeType.WALL):
                Vor[i,j] = 2 * ( Psi[i,j] - Psi[i-1,j] )/ dx2

            #bottom wall
            if ( j == 0 and nodes_type[i,j] == NodeType.WALL ):
                Vor[i,j] = 2 *( Psi[i,j] - Psi[i,j+1] )/ dy2

            #top wall
            if ( nodes_type[i,j] == NodeType.LID ):
                Vor[i,j] = 2 *( Psi[i,j] - Psi[i,j-1] )/ dy2  -  2 * u0 / dx


'''advances vorticity equation using RK4'''
def advanceRK4():

    #compute the four terms of RK4
    RK = RHS( VorOld )
    Vor1 = VorOld + 0.5*dt*RK

    RK1 = RHS( Vor1 )
    Vor2 = VorOld + 0.5*dt*RK1

    RK2 = RHS( Vor2 )
    Vor3 = VorOld + dt*RK2

    RK3 = RHS( Vor3 )
    VorNew = VorOld + (dt/6.0)* (RK + 2*RK1 + 2*RK2 + RK3) 
    
    return VorNew


'''plot distribution of nodes'''
def plotNodes():

    for i in range( 0, MB ):
        for j in range( 0, NB ):

            if nodes_type[i, j] == NodeType.WALL:
                color ='k'
            elif nodes_type[i, j] == NodeType.LID:
                 color ='m'
            else:
                color ='g'
            plt.scatter(pos_x[i, j], pos_y[i, j],c=color)

    plt.show()

'''PrintOutput'''     
def printDate():

    a = np.reshape(pos_x,MB*NB)
    b = np.reshape(pos_y,MB*NB)
    c =  np.reshape(u,MB*NB)
    d =  np.reshape(v,MB*NB)
    e =  np.reshape(Vor,MB*NB)
    f =  np.reshape(Psi,MB*NB)
    g =  np.reshape(nodes_type,MB*NB)  
    
    np.savetxt('.\OutPut\Data'+str("%.5f" % (tStep*dt))+'.txt', np.column_stack((a, b, c, d, e, f, g)),
               fmt="%2.3f", delimiter=" , ")

'''main program'''                     
# parameters of geometry
Length = 1     # length of computational domain
Height = 1      # height of computational domain
MB = 41         # number of the nodes in x direction
NB = 41         # number of the nodes in y direction

dx = Length / (MB-1)
dy = Height / (NB-1)

# parameters of fluid
nu = 0.1      # kinematic viscosity 
u0 = 1          # inlet velocity

# parameters for computation
time_total = 10  # total time
dt =0.001        # time step size
recordN =  100      # time steps for recording data 
maxiter = 500   # maximum number of iteration

# setting folder for saving data
dir_name = "OutPut"

if not os.path.exists(dir_name):    #先確認資料夾是否存在
    os.makedirs(dir_name)
else:
    shutil.rmtree(dir_name)  
    os.makedirs(dir_name)



# generate geometry
pos_x, pos_y, nodes_type = makeGeometry()

# setting initial condition of Psi
Psi, Vor, u, v = IntPsi()




# start iteration
print ("Starting main loop")

for tStep in range ( int(time_total/dt) ):
    
    # Update Solutions
    VorOld = Vor
    PsiOld = Psi

    for it in range(maxiter+1):
        tempPsi = Psi
        tempVor = Vor

        # solve Psi
        computePsi()

        # update u and v
        computeVel()

        # advance w
        Vor = advanceRK4()
        VorBou()

        #check for convergence
        error = np.linalg.norm(Vor-tempVor) / np.linalg.norm(Vor) 
        if (error<=1e-5):
            print('computation is convergent! ')
            print('convergence in iter : ',it,' , ', error )
            break

        elif (it >= maxiter):
            print ('Warning: solutions of Psi are not convegent')
            print('convergence in iter : ',it,' , ', error )
            print ('Program end...')
            exit()

        else:
            print('convergence in iter : ',it,' , ', error )

    if (tStep % recordN == 0):
        printDate()
    
    if ((tStep+1)*dt == time_total):
        tStep += 1
        printDate()

    fig, ax = plt.subplots(figsize=(10, 8), dpi=100, facecolor='w', edgecolor='k')



# plotNodes()