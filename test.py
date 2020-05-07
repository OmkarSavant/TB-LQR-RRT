import numpy as np 
import networkx as nx
from rtree import index
# import rospy

class BallisticMotion():
    def __init__(self):
        """
        Class to handle the dynamics of the ball
        """

    def collision():
        """
        Updates dynamics after a collision occurs
        """
        return

    def get_position(self, t):
        """
        Returns the position of the ball given its dynamics
        """
        return np.array([0, t, 1.0, 0, 0, 0])


class TBRRT():
    def __init__(self, start, goal, mode = 'manipulator'):
        """
        :param start: robot arm starting position
        :param goal: goal position
        """
        
        # Setup RTree
        dim = len(start)
        p = index.Property()
        p.dimension = dim
        self.Tree = index.Index(interleaved=True, properties=p)

        # Setup Graph
        self.G = nx.Digraph()

        # Store variables
        self.start = start
        self.goal = goal #TODO: How to select the goal such that it interfaces with the ball flight?
        self.mode = mode

        if (self.mode == 'manipulator'):
            self.min_state = np.array([0, -5*np.pi/6, -5*np.pi/6, -20, -20, -20])
            self.max_state = np.array([2*np.pi, 5*np.pi/6, 5*np.pi/6, 20, 20, 20])



    def sample(self):
        """
        Sample the joint space
        """
                    
        return np.random.uniform(self.min_state, self.max_state)

    def nearest_neighbor(self, xRand):
        """
        Find the node in the graph which requires the smallest magnitude of u to get to from the random state.
        TODO: Given the control limits, can we calculate state bounds and use this as a heuristic? Or just check within a radius. Or
        terminate if control is less than some value?
        Approximate nearest neighbors
        - Can also use sqrt(N) samples
        """

        # DON'T LOOK AT TIME FOR NEAREST NEIGHBOR, BUT RETRUN FULL STATE-TIME VECTOR!

        min_state = ()
        min_control = np.inf

        for state in self.G.nodes():
            
            #compute the u vector required to reach this state

            #f(xk,uk) = M_inverse * (u - C(x1,x2)*x2 - G(x1))
            #xK+1 = xk + f(xk, uk)*dt 

            J1 = np.array([[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]])
            
            
            f = (state - xRand)/self.dt

            u = M * f + self.C(xRand)*x2 + self.G(xRand)
            
            if control < min_control:
                
                min_state = state
                min_control = control
            
        return min_state #This should be the state-time vector

    def C(self, state):
        """
        """
        return

    def steer(self, xi, xRand):
        """
        Select ui to move from xi toward xrand
        """
        return ui

    def dynamics(self, xi, ui):
        """
        :param xi: xi is the full state time vector, but only the state will be used
        returns state time_vector after applying dynamics as described by f(xi, ui)
        """
        #TODO
        return velocities?

    def is_reachable(self, xi):
        """
        Checks if xi is reachable, i.e. not self collisions or joint limits
        """
        return True

    def search(self, max_samples):
        """
        Perform TB-RRT Algorithm
        :param max_samples: Number of samples until termination 
        """
        for k in range(max_samples):

            xRand = self.sample()
            xi = self.nearest_neighbor(xRand)
            xRand = self.sample()
            ui = self.steer(xi, xRand)
            xi_1 = xi + np.array([self.dynamics(xi, ui), 1]) * self.dt

            if self.is_reachable(xi_1):
                

            
            
            
            
        # for k = 1 --> K: 
            #xRand = random state-time vector
                #xi <-- nearest_neighbor(T,xrand)
                #select ui to move from xi towards xRand
                #x_i+1 <-- xi + [f(xi,ui) 1]^T * delta_T
                
                #if xi+1 is in set Xfree:

                    #add xi+1 to Tree 
                    # add edge (xi,xi+1) to T  with weight of control ui

                    
                    #if metric function rho(xi+1, X_goal) < eps:
                        
                        #calculate path P in Tree from xo to xi+1

                        #return P


"""

- "Nearest neighbor" is the state in the tree that requires the least amount of control effort to get to x_rand. 
Questions: How do determine the control cost between two states? For example if x_i is moving in the opposite direction of x_rand, then what
is the cost to make it happen? If we assume instantaneous acceleration. 
How do we reduce the the computational complexity? Pruning? Distance or velocity direction heuristics?
We probably don't need Rtree...

- When random samples are found, they are done in the position, velocity space. When they are stored in the tree, 
   they are tagged with the appropriate time step by setting it to the parent's time + delta_T

- Edge costs are the norm of controls applied when transitioning between two states. This edge should map to some 
dictionary which contains the 6d vector of actual joint torques during this transitions 

- only using their metric function for assessing proximity to goal 

"""




if __name__=='__main__':
    
    #state = theta (degrees), theta_dot, time
    
    #Start position 
    start = np.array([0,0])

    #end position
    end = np.array([15,0.5])

    #intialize Tree with Xo (starting point)
    rrt = TBRRT(start,end)
    path = rrt.search(10000)            