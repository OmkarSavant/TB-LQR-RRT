import numpy as np
import pdb
import matplotlib.pyplot as plt
import scipy.linalg as sp_linalg

class Node:
    def __init__(self, state_time = np.zeros(3, dtype=float)):
        
        self.state_time = state_time
        #lqr params
        self.u = 0
        self.parent = None
        self.children = []

class TBRRT():
    def __init__(self, arm, start, goal, dt = 0.01):
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

        self.maxU = 500
        self.minUs = -self.maxUs

    def sample(self):
        """
        Sample the joint space
        """

        #compute the relevant range between the starting position and ending position that should be sampled 
        rand_state = np.random.uniform(self.min_state, self.max_state)
        
        return rand_state

    def compute_dist(self, s1, s2):
        """
        Return the norm between 2 state-time vectors. Note to wrap the angle difference
        """
        return sp_linalg.norm(self.wrap(s1, s2))
        

    def nearest_neighbor(self, xRand):
        """
        Find the node in the graph which requires the smallest magnitude of u to get to from the random state.
        """
        # TODO: Make this more efficient?
        #within a neighborhood of XRand, determine the lowest cost to go
        minCost = np.inf
        minNode = None

        for node in self.Tree:

            cost = self.compute_dist(node.state_time[0:6], xRand)

            if cost < minCost:
                minNode = node
                minCost = cost

        return minNode

    def steer(self, xi, xRand, k = 10):
        """
        Select ui to move from xi toward xrand
        """

        #test k different controls 
        minDist = np.inf
        bestControl = np.zeros(3)
        xi_1 = None 

        for _ in range(k):
            
            #smarter sampling in the right direction? 
            u_rand = np.random.uniform(self.minUs, self.maxUs)
            nextState = self.simulate_system(xi.state_time[0:6], u_rand)
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
        #TODO: Check if any of the links intersect...
        return True

    def add_edge(self, parent, child, ui):
        parent.children.append(child)

        child.parent = parent
        
        self.edge_costs[(parent, child)] = ui

    def near_goal(self, xi_1):
        
        dist_norm = np.linalg.norm(xi_1[0:3] - self.goal[0:3])
        vel_norm = np.linalg.norm(xi_1[3:7] - self.goal[3:7])
        time_norm = np.linalg.norm(xi_1[-1] - self.goal[-1])

        if (dist_norm < self.eps[0]) and (vel_norm < self.eps[1]) and (time_norm < self.eps[2]):
            return True
        else:
            return False

    def get_path(self, end_node):
        
        path = []
        node = end_node
        
        while node.parent != None:
            path.append(node.parent.state_time)
            node = node.parent

        path.append(node)
        return path.reverse()

    def wrap(self, target, source):
        """
        Returns wrapped angle
        """
        
        #if you are only checking angular difference
        if len(target) == 3:
            return np.arctan2(np.sin(target[0:3]-source[0:3]), np.cos(target[0:3]-source[0:3]))
        else:
            return np.append(np.arctan2(np.sin(target[0:3]-source[0:3]), np.cos(target[0:3]-source[0:3])), target[3:] - source[3:])
    
    def simulate_system(self, state, input, time_step = 0.01): 
        """
        Simulate system. Only works on state, NOT state time
        """

        #set arm position to x
        self.arm.reset(q=state[0:3],dq=state[3:6])

        #apply the control signal
        self.arm.apply_torque(input,time_step)

        #get the next step from the arm 
        xnext = np.append(np.copy(self.arm.q),np.copy(self.arm.dq))

        return xnext

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
                    path = self.get_path(xi_1_node)
                    return path

                if len(self.Tree) % 1000 == 0:
                    #TODO: This doesnt consider time...
                    nearest_to_goal = self.nearest_neighbor(self.goal[0:6])
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
    
    subfolder = 'three_link'
    arm_name = 'arms.%s.%s' % (subfolder, 'arm')
    arm_module = importlib.import_module(name=arm_name)
    arm = arm_module.Arm(dt=0.02)

    #state = theta (degrees), theta_dot, time
    
    #Start position 
    start = np.array([0., 0., 0., 0., 0., 0., 0.])

    #end position
    end = np.array([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1]) # This is random...

    #intialize Tree with Xo (starting point)
    rrt = TBRRT(arm, start,end)

    path = rrt.search(100000)