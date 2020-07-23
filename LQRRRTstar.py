import numpy as np
import pdb
import matplotlib.pyplot as plt
import scipy.linalg as sp_linalg
import control
from itertools import islice
import time
from mpl_toolkits import mplot3d

import math

class Node:
    def __init__(self, state_time = np.zeros(3, dtype=float)):
        
        self.state_time = state_time
        #lqr params
        self.u = 0
        self.parent = None
        self.children = []
        self.cost_so_far = 0

class MCRRT():
    def __init__(self, start, goal, dt = 0.01):
        
        """
        :param start: robot arm starting position
        :param goal: goal position
        """

        #for k nearest
        self.k = 10 

        self.Tree = []
        self.edge_costs = {} # (parent, child) -> cost

        self.start = start
        self.goal = goal

        # Store variables
        self.dt = dt

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
        

    def nearest_neighbor(self, xRand, k = 1):
        """
        Find the node in the graph which requires the smallest magnitude of u to get to from the random state.
        """
        # TODO: Make this more efficient?
        #within a neighborhood of XRand, determine the lowest cost to go

        nodeDict = {}

        for node in self.Tree:

            if node.state_time[-1] > self.goal[-1]:
                continue

            #xRand = xRand.reshape(2, 1)
            v_minus_x = np.matrix(node.state_time[0:-1] - xRand[0:-1])
            cost = np.matmul(np.matmul(v_minus_x, self.S), v_minus_x.T)
            nodeDict[node] = cost

        #sort this dictionary and output the 10 lowest values 
        orderedNodes = {k: v for k, v in sorted(nodeDict.items(), key=lambda item: item[1])}
        
        k_orderedNodes = [items[0] for items in list(islice(orderedNodes.items(), k))]

        return k_orderedNodes
        
    
    def TB_nearest_neighbor(self, xRand, k, timeStep):
        """
        Find the node in the graph which requires the smallest magnitude of u to get to from the random state.
        """
        # TODO: Make this more efficient?
        #within a neighborhood of XRand, determine the lowest cost to go

        nodeDict = {}

        #print('new tb call:' + str(xRand[-1]) + 'timestep: ' + str(timeStep))
        for node in self.Tree:

            #print(node.state_time[-1])

            #if the time step is not from the previous time step, not interested
            if math.isclose(node.state_time[-1],(xRand[-1] + timeStep*self.dt)) == False:
                continue
            
            #print('new child' + str(node.state_time[-1]))
            #xRand = xRand.reshape(2, 1)
            v_minus_x = np.matrix(node.state_time[0:-1] - xRand[0:-1])
            cost = np.matmul(np.matmul(v_minus_x, self.S), v_minus_x.T)
            nodeDict[node] = cost

        #sort this dictionary and output the 10 lowest values 
        orderedNodes = {k: v for k, v in sorted(nodeDict.items(), key=lambda item: item[1])}
        
        k_orderedNodes = [items[0] for items in list(islice(orderedNodes.items(), k))]

        if len(k_orderedNodes) == 0 and timeStep == -1:
            pdb.set_trace()

        return k_orderedNodes

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

    def calc_lqr_cost(self,x1,x2):

        #determine LQR cost of going from x1 to x2
        
        x_del = x1.state_time[0:-1] - x2.state_time[0:-1]

        u = np.array([-self.K * x_del.reshape(2, 1)]).squeeze()    
       
        return u
        

    def calc_lqr(self):
         K, S, _ = control.lqr(self.A, self.B, self.Q, self.R)

         return K, S

    def chooseParent(self, X_near, xRand):
        """
        X_near: list of k nearest NODES
        xRand: randoom STATE vector (not a node)
        """
 
        minCost = np.inf
        x_min = None
        sigma_min = None

        for x_near in X_near:

            sigma, u = self.steer(x_near,xRand)

            if x_near.cost + u < minCost:
                minCost = x_near.cost + u
                x_min = x_near
                sigma_min = sigma

        return x_min, sigma_min
    
    def plot_tree(self, path):
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
    
        #theta
        xdata = [node.state_time[0] for node in self.Tree]

        #theta_dot
        ydata = [node.state_time[1] for node in self.Tree]

        #time
        zdata = [node.state_time[2] for node in self.Tree]

        ax.scatter3D(xdata,ydata,zdata, color = 'blue',alpha = 0.3)
        ax.set_xlabel('theta (rad)')
        ax.set_ylabel('theta_dot (rad/s)')
        ax.set_zlabel('time (s)')
    
        #plot nodes on path
        x_path = [node.state_time[0] for node in path]
        y_path = [node.state_time[1] for node in path]
        z_path = [node.state_time[2] for node in path]

        ax.scatter3D(x_path,y_path,z_path, color = 'red', alpha=1.0)
        
        # Plot goal hyper plane....


        # Plot start node

        plt.savefig('tree'+str(int(time.time())))

    def search(self, max_samples = 10000):
        """
        Perform TB-RRT Algorithm
        :param max_samples: Number of samples until termination 
        """

        while (len(self.Tree) < max_samples):
            
            xRand = self.sample()
            xNearest = self.nearest_neighbor(xRand, k = 1)[0]
            xNew, cost = self.steer(xNearest, xRand)
            
            # Prepare xNew as node
            xNew_stateTime = np.append(xNew, xNearest.state_time[-1] + self.dt)
            xNew = Node(xNew_stateTime)
            xNew.cost_so_far = xNearest.cost_so_far + cost # Check this

            #should be collissionFree 
            if self.is_reachable(xNew):
                k = np.minimum(np.maximum(int(2*np.e*np.log(len(self.Tree))), 1), len(self.Tree))
                
                #list of nodes
                # Choose best parent
                # from the past
                xNear = self.TB_nearest_neighbor(xNew.state_time, k, -1.0)

                xMin = None
                cMin = np.inf

                for xn in xNear:

                    tmp_cost = xn.cost_so_far + self.calc_lqr_cost(xn, xNew)

                    if tmp_cost < cMin:
                        xMin = xn
                        cMin = tmp_cost

                # Add to tree
                xNew.cost_so_far = cMin
                self.Tree.append(xNew)
                
                # Add edge
                xNew.parent = xMin
                xMin.children.append(xNew)
                xNew.state_time[-1] = xNew.parent.state_time[-1] + self.dt
                
                # Rewire
                xNear = self.TB_nearest_neighbor(xNew.state_time, k, 1.0)

                #xNear should be from the future
                for xn in xNear:
                    
                    #another get_cost call 
                    tmp_cost = xNew.cost_so_far + self.calc_lqr_cost(xNew, xn)
                    if tmp_cost < xn.cost_so_far:

                        xn.parent.children.remove(xn)
                        xn.parent = xNew    
                        xNew.children.append(xn)
                        xn.cost_so_far = tmp_cost
                        xNew.state_time[-1] = xNew.parent.state_time[-1] + self.dt

                if len(self.Tree) % 20 == 0:
                    print(len(self.Tree))

        # Find path after some number of iterations
        for node in self.Tree:
            if self.near_goal(node.state_time):
                path, total_control = self.get_path(node)
                self.plot_tree(path)
                return path, total_control
            
            #even if there's no path, get the closest node to the end 
            else:
                pass

        
        print("No path found")
        pdb.set_trace()

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

    path = rrt.search(1000) #yeah
    pdb.set_trace()