import numpy as np
import pdb
import matplotlib.pyplot as plt
import scipy.linalg as sp_linalg
import control

class Node:
    def __init__(self, state_time = np.zeros(3, dtype=float)):
        
        self.state_time = state_time
        #lqr params
        self.u = 0
        self.parent = None
        self.children = []

class MCRRT():
    def __init__(self, start, goal, dt = 0.01):
        
        """
        :param start: robot arm starting position
        :param goal: goal position
        """
        self.Tree = []
        self.edge_costs = {} # (parent, child) -> cost

        self.start = start
        self.goal = goal

        # Store variables
        self.dt = dt

        self.goal = goal
        self.eps = np.array([np.deg2rad(2.8), 0.2, 0.03])

        # Add start node to the tree
        start_node = Node(start)
        self.Tree.append(start_node)

        self.min_state = np.array([0,-6])
        self.max_state = np.array([np.pi/2,6])

        self.maxU = 2
        self.minU = -self.maxU

        self.m = 2
        self.l = 0.2

        self.A = np.array([[0,1],[0,0]])
        self.B = np.array([[0],[3/(self.m * self.l**2)]])
        self.Q = np.array([[10,0],[0,10]])
        self.R = np.array([1])

        self.K, self.S = self.calc_lqr()
        
        self.create_smart_sample()
        
    def create_smart_sample(self):
        
        #determine distance between start and goal
        
        smart_min_theta = np.minimum(np.minimum(self.goal[0],self.start[0]) - 0.1, self.min_state[0])
    
        smart_max_theta = np.maximum(np.maximum(self.goal[0],self.start[0]) + 0.1,self.max_state[0])

        self.smart_min_state = np.array([smart_min_theta,self.min_state[1]])
        self.smart_max_state = np.array([smart_max_theta,self.max_state[1]])
        
    def sample(self):
        """
        Sample the joint space
        """

        #compute the relevant range between the starting position and ending position that should be sampled 
        return np.random.uniform(self.smart_min_state, self.smart_max_state)

    def compute_dist(self, s1, s2):
        """
        Return the norm between 2 state-time vectors. Note to wrap the angle difference
        """
        return sp_linalg.norm(s1 - s2)
        

    def nearest_neighbor(self, xRand):
        """
        Find the node in the graph which requires the smallest magnitude of u to get to from the random state.
        """
        # TODO: Make this more efficient?
        #within a neighborhood of XRand, determine the lowest cost to go
        minCost = np.inf
        minNode = None

        for node in self.Tree:

            #xRand = xRand.reshape(2, 1)
            v_minus_x = np.matrix(node.state_time[0:-1] - xRand)

            cost = np.matmul(np.matmul(v_minus_x, self.S), v_minus_x.T)

            if cost < minCost:
                minNode = node
                minCost = cost

        return minNode

    def steer(self, xi, xRand):
        """
        Select ui to move from xi toward xrand
        """

        x_del = xi.state_time[0:-1] - xRand
        u = np.array([-self.K * x_del.reshape(2, 1)]).squeeze()    
        
        #print(u)

        if u > 2: 
            u = 2
        if u < -2:
            u = -2

        nextState = self.simulate_system(xi.state_time[0:-1], u)
        
        return nextState, u

    def is_reachable(self, xi):
        """
        Checks if xi is reachable, i.e. not self collisions or joint limits
        """
        return True

    def add_edge(self, parent, child, ui):
        parent.children.append(child)
        child.parent = parent
        self.edge_costs[(parent, child)] = ui

    def near_goal(self, xi_1):
        
        dist_norm = np.linalg.norm(xi_1[0:1] - self.goal[0:1])
        vel_norm = np.linalg.norm(xi_1[1:2] - self.goal[1:2])
        time_norm = np.linalg.norm(xi_1[-1] - self.goal[-1])

        if (dist_norm < self.eps[0]) and (vel_norm < self.eps[1]) and (time_norm < self.eps[2]):
            return True
        else:
            return False

    def get_path(self, end_node):
        
        path = []
        node = end_node
        control_cost = 0
        
        while node.parent != None:
            path.append(node.parent.state_time)
            node = node.parent
            control_cost += abs(node.u)

        path.append(node)
        return path.reverse(), control_cost
    
    def simulate_system(self, state, input, time_step = 0.01): 
        """
        Simulate system. Only works on state, NOT state time
        """

        #get the next step from the arm 
        f_x_u = np.array([state[1], 3*input/(self.m*self.l**2)])
        xnext = state + f_x_u * time_step

        return xnext

    def calc_lqr(self):
         K, S, _ = control.lqr(self.A, self.B, self.Q, self.R)

         return K, S

    def search(self, max_samples = 10000):
        """
        Perform TB-RRT Algorithm
        :param max_samples: Number of samples until termination 
        """

        while (len(self.Tree) < max_samples):
            
            xRand = self.sample()

            xi = self.nearest_neighbor(xRand)

            # xi = self.Tree[np.random.randint(len(self.Tree))]

            xi_1, xi_u = self.steer(xi, xRand)

            if self.is_reachable(xi_1):

                xi_1_stateTime = np.append(xi_1, xi.state_time[-1] + self.dt)
                
                if xi_1_stateTime[-1] > self.goal[-1]:
                    continue

                xi_1_node = Node(xi_1_stateTime)
                
                xi_1_node.u = xi_u
                               
                self.Tree.append(xi_1_node)

                if len(self.Tree) % 50 == 0:
                    print(len(self.Tree))
                self.add_edge(xi, xi_1_node, np.linalg.norm(xi_u))
                
                if self.near_goal(xi_1_stateTime):
                    path, control_cost = self.get_path(xi_1_node)
                    print("path found")
                    print(control_cost)
                    return path, control_cost

                # Debugging
                if len(self.Tree) % 5000 == 0:
                    #TODO: This doesnt consider time...
                    nearest_to_goal = self.nearest_neighbor(self.goal[0:-1])
                    longest_path = -np.inf
                    for node in self.Tree:
                        t = node.state_time[-1]
                        if t > longest_path:
                            longest_path = t
                    print(longest_path)
                    pdb.set_trace()
        
        print("No path found")

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
    
    #Start position 
    start = np.array([np.pi/4., 0., 0.])

    #end position
    end = np.array([np.deg2rad(10), 2, 0.3]) 

    #intialize Tree with Xo (starting point)
    rrt = MCRRT(start, end)

    path = rrt.search(20000)