import numpy as np
import pdb
import matplotlib.pyplot as plt
import scipy.linalg as sp_linalg
import control

#may 10 

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

    def sample(self):
        """
        Sample the joint space
        """

        #compute the relevant range between the starting position and ending position that should be sampled 
        return np.random.uniform(self.min_state, self.max_state)

    def compute_dist(self, s1, s2):
        """
        Return the norm between 2 state-time vectors. Note to wrap the angle difference
        """
        return sp_linalg.norm((s1 - s2) * np.array([4 , 1]))
        

    def nearest_neighbor(self, xRand):
        """
        Find the node in the graph which requires the smallest magnitude of u to get to from the random state.
        """
        # TODO: Make this more efficient?
        #within a neighborhood of XRand, determine the lowest cost to go
        minCost = np.inf
        minNode = None

        for node in self.Tree:

            cost = self.compute_dist(node.state_time[0:-1], xRand)

            if cost < minCost:
                minNode = node
                minCost = cost

        return minNode

    def steer(self, xi, xRand, k = 20):
        """
        Select ui to move from xi toward xrand
        """

        #test k different controls 
        minDist = np.inf
        bestControl = 0
        xi_1 = None 

        for _ in range(k):
            
            #smarter sampling in the right direction? 
            u_rand = np.random.uniform(self.minU, self.maxU)
            nextState = self.simulate_system(xi.state_time[0:-1], u_rand)
            next_state_dist = self.compute_dist(nextState, xRand)

            if next_state_dist < minDist:
                
                minDist = next_state_dist
                xi_1 = nextState
                bestControl = u_rand

        return xi_1, minDist, bestControl

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
        path.append(end_node)
        
        while node.parent != None:
            path.append(node.parent)
            node = node.parent
            control_cost += abs(node.u)

        path.append(node)
        path.reverse()

        return path, control_cost
    
    def simulate_system(self, state, input, time_step = 0.01): 
        """
        Simulate system. Only works on state, NOT state time
        """

        #get the next step from the arm 
        f_x_u = np.array([state[1], 3*input/(self.m*self.l**2)])

        xnext = state + f_x_u * time_step

        return xnext

    # def calc_lqr(self):
    #     self.K, self.S, _ = control.lqr(self.A, self.B, self.Q, self.R)

    def search(self, max_samples = 10000):
        """
        Perform TB-RRT Algorithm
        :param max_samples: Number of samples until termination 
        """
        while (len(self.Tree) < max_samples):
            
            xRand = self.sample()

            xi = self.nearest_neighbor(xRand)

            # xi = self.Tree[np.random.randint(len(self.Tree))]

            xi_1, minDist, bestControl = self.steer(xi, xRand)

            if self.is_reachable(xi_1):

                xi_1_stateTime = np.append(xi_1, xi.state_time[-1] + self.dt)
                
                if xi_1_stateTime[-1] > self.goal[-1]:
                    continue

                xi_1_node = Node(xi_1_stateTime)
                
                xi_1_node.u = bestControl
                               
                self.Tree.append(xi_1_node)

                if len(self.Tree) % 10 == 0:
                    print(len(self.Tree))
                self.add_edge(xi, xi_1_node, np.linalg.norm(bestControl))
                
                if self.near_goal(xi_1_stateTime):
                    path, control_cost = self.get_path(xi_1_node)
                    print(control_cost)
                    return path, control_cost, len(self.Tree)

                # Debugging
                if len(self.Tree) % 4200 == 0:
                    #TODO: This doesnt consider time...
                    return None, None, len(self.Tree)
#                    nearest_to_goal = self.nearest_neighbor(self.goal[0:-1])
#                    longest_path = -np.inf
#                    for node in self.Tree:
#                        t = node.state_time[-1]
#                        if t > longest_path:
#                            longest_path = t
#                    print(longest_path)
        
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

- only using their metric functi
on for assessing proximity to goal 

"""

if __name__=='__main__':
    
#    #Start position 
#    start = np.array([np.pi/2., 0., 0.])
#
#    #end position
#    end = np.array([np.deg2rad(10), 2, 0.5]) 
#
#    #intialize Tree with Xo (starting point)
#    rrt = MCRRT(start, end)
#
#    path, cost = rrt.search(20000)
    
    costs = []
    its = []
    
    #Start position 
    start = np.array([np.pi/2., 0., 0.])
        
    #end position
    end = np.array([np.deg2rad(10), 2, 0.5]) 
    
    for i in range(10):
    
        #intialize Tree with Xo (starting point)
        rrt = MCRRT(start, end)
        
        path, cost, it = rrt.search(20000)
        costs.append(cost)
        its.append(it)