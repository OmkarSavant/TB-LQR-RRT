import numpy as np
import control
import pdb
import matplotlib.pyplot as plt
import importlib

import scipy.linalg as sp_linalg

# class BallisticMotion():
#     def __init__(self):
#         """
#         Class to handle the dynamics of the ball
#         """

#     def collision(self):
#         """
#         Updates dynamics after a collision occurs
#         """
#         return

#     def get_position(self, t):
#         """
#         Returns the position of the ball given its dynamics
#         """
#         return np.array([0, t, 1.0, 0, 0, 0])

class Node:
    def __init__(self, state_time = np.zeros(7, dtype=float)):
        
        self.hasLQR = False
        self.state_time = state_time

        #lqr params
        self.K = None
        self.S = None
        self.u = np.zeros(3)

        self.parent = np.nan
        self.child = np.nan

class TBRRT():
    def __init__(self, arm, start, goal):
        """
        :param start: robot arm starting position
        :param goal: goal position
        """
        self.arm = arm
        self.Tree = []
        self.edge_costs = {} # (parent, child) -> cost

        # Store variables
        self.dt = arm.dt

        self.goal = goal
        self.eps = np.array([2.8, 0.2, 0.03])

        self.Q = np.zeros((6,6))
        self.Q[:3,:3] = np.eye(3) * 1000
        self.R = np.eye(3) * 0.001

        # Add start node to the tree
        start_node = Node(start)
        self.gen_lqr_new(start_node)
        self.Tree.append(start_node)

        self.min_state_time = np.array([0,0,0,-1,-1,-1,0])
        self.max_state_time = np.array([2*np.pi,2*np.pi,2*np.pi,1,1,1,self.goal[-1]])

        self.minUs = np.array([-1000,-1000,-1000])
        self.maxUs = np.array([1000,1000,1000]) * 100

    def sample(self):
        """
        Sample the joint space
        """
        rand_state_time = np.random.uniform(self.min_state_time, self.max_state_time)
        xRand = Node(rand_state_time)

        self.gen_lqr_new(xRand)
        
        return xRand

    def nearest_neighbor(self, xRand):
        """
        Find the node in the graph which requires the smallest magnitude of u to get to from the random state.
        """

        #within a neighborhood of XRand, determine the lowest cost to go
        minCost = np.inf
        minNode = None

        for node in self.Tree:

            #filter for nodes that are within a 15 degree ball for all three angular states to ensure linearlity of LQR approx
            #handle wrapping and absolute, convert from rad to degree
            # https://stackoverflow.com/questions/1878907/the-smallest-difference-between-2-angles
            # if abs(np.arctan2(np.sin(node.state_time[0] - xRand.state_time[0]), np.cos(node.state_time[0] - xRand.state_time[0]))) > 0.262:
            #     continue

            # if abs(np.arctan2(np.sin(node.state_time[1] - xRand.state_time[1]), np.cos(node.state_time[1] - xRand.state_time[1]))) > 0.262:
            #     continue

            # if abs(np.arctan2(np.sin(node.state_time[2] - xRand.state_time[2]), np.cos(node.state_time[2] - xRand.state_time[2]))) > 0.262:
            #     continue

            #check that all 3 angles are w/in 15 degrees of the node in consideration before weighing its LQR cost
            
            if np.any(self.wrap(node.state_time[0:3], xRand.state_time[0:3]) > np.ones(3)*np.deg2rad(15.0)):
                #find the cost to go from the neighbor to xRand
                v_minus_x = np.matrix(self.wrap(node.state_time[0:-1], xRand.state_time[0:-1]))

                cost = np.matmul(np.matmul(v_minus_x, xRand.S), v_minus_x.T)

                if cost < minCost:
                    minNode = node
                    minCost = cost

        return minNode #This should be the state-time vector

    def steer(self, xi, xRand):
        """
        Select ui to move from xi toward xrand
        """
        #use the K matrix from xi, and the difference in state should be (xRand-xi)
        x_del = self.wrap(xi.state_time[0:-1] , xRand.state_time[0:-1]).reshape(6,1)

        ui = np.array(-1 * np.matmul(xi.K, x_del)).reshape((3,))

       # ui_clipped = np.clip(ui,self.minUs,self.maxUs)

        #apply control limits to ui 
       # ui_clipped = np.clip(ui,self.minUs, self.maxUs)

        xi_1 = self.simulate_system(xi, xi.state_time[0:6], ui, self.dt)

        return xi_1, ui

    def is_reachable(self, xi):
        """
        Checks if xi is reachable, i.e. not self collisions or joint limits
        """
        #TODO: Check if any of the links intersect...
        return True

    def add_edge(self, parent, child, ui):
        parent.child = child
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

    def gen_lqr_new(self, node):
        
        node.hasLQR = True

        u = np.zeros((3,))

        x = node.state_time[0:-1]

        eps = 0.00001  # finite difference epsilon
        #----------- compute xdot_x and xdot_u using finite differences --------
        # NOTE: here each different run is in its own column
        x1 = np.tile(x, (6,1)).T + np.eye(6) * eps
        x2 = np.tile(x, (6,1)).T - np.eye(6) * eps
        uu = np.tile(u, (6,1))

        f1 = self.plant_dynamics(x1, uu)
        f2 = self.plant_dynamics(x2, uu)

        xdot_x = (f1 - f2) / 2 / eps
   
        xx = np.tile(x, (3,1)).T 
        u1 = np.tile(u, (3,1)) + np.eye(3) * eps
        u2 = np.tile(u, (3,1)) - np.eye(3) * eps

        f1 = self.plant_dynamics(xx, u1)
        f2 = self.plant_dynamics(xx, u2)

        xdot_u = (f1 - f2) / 2 / eps

        self.A = xdot_x
        self.B = xdot_u

        node.S = sp_linalg.solve_discrete_are(self.A, self.B, self.Q, self.R)
        node.K = np.dot(np.linalg.pinv(self.R + np.dot(self.B.T, np.dot(node.S, self.B))), np.dot(self.B.T, np.dot(node.S, self.A)))

        #node.K, node.S, _ = control.lqr(self.A, self.B, self.Q, self.R)

    def gen_lqr(self, node):

        """ 
        Create a local estimate of the K and S matrices

        Linearize about x, with u = 0 (as per paper)
        """
        node.hasLQR = True

        control_signal = node.u
        currState = node.state_time[0:6]

        eps = 0.00001

        time_step = self.dt
        
        A = np.zeros([len(currState),len(currState)])
        for i in range(len(currState)):
            x = currState.copy()
            x[i] += eps
            x_inc = self.simulate_system(node, x, control_signal, time_step)
            x = currState.copy()
            x[i] -= eps
            x_dec = self.simulate_system(node, x, control_signal, time_step)
            A[:, i] = self.wrap(x_inc, x_dec) / (2*eps)

        B = np.zeros([len(currState), len(control_signal)])
        for i in range(len(control_signal)):
            u = control_signal.copy()
            u[i] += eps
            x_inc = self.simulate_system(node, currState, u, time_step)
            u = control_signal.copy()
            u[i] -= eps
            x_dec = self.simulate_system(node, currState, u, time_step)
            B[:,i] = self.wrap(x_inc, x_dec) / (2*eps)
        
        node.K, node.S, _ = control.lqr(A, B, self.Q, self.R)
    
    def wrap(self, target, source):
        """
        Returns wrapped angle
        """
        
        #if you are only checking angular difference
        if len(target) == 3:
            return np.arctan2(np.sin(target[0:3]-source[0:3]), np.cos(target[0:3]-source[0:3]))
        else:
            return np.append(np.arctan2(np.sin(target[0:3]-source[0:3]), np.cos(target[0:3]-source[0:3])), target[3:] - source[3:])
    
    def simulate_system(self, node, state, input, time_step = 0.01): 
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
    
    def plant_dynamics(self, x, u):

        """ Simulate the arm dynamics locally. """

        xnext = np.zeros((x.shape))

        for ii in range(x.shape[1]):
            # set the arm position to x
            self.arm.reset(q=x[:3, ii], 
                          dq=x[3:6, ii])

            # apply the control signal
            # TODO: should we be using a constant timestep here instead of arm.dt?
            # to even things out when running at different dts? 

            self.arm.apply_torque(u[ii], self.arm.dt)
            # get the system state from the arm
            xnext[:,ii] = np.hstack([np.copy(self.arm.q), 
                                   np.copy(self.arm.dq)])

        return xnext

    def search(self, max_samples = 10000):
        """
        Perform TB-RRT Algorithm
        :param max_samples: Number of samples until termination 
        """
        while (len(self.Tree) < max_samples):
            xRand = self.sample()
            xi = self.nearest_neighbor(xRand)
            
            if(xi == None):
                # print("Too far")
                continue

            xi_1, ui = self.steer(xi, xRand)

            if (np.any(abs(ui) > self.maxUs)):
                #print("Too much control")
                continue

            if self.is_reachable(xi_1):
                xi_1 = np.append(xi_1, xi.state_time[-1]+self.dt)
                xi_1_node = Node(xi_1)
                xi_1_node.u = ui

                self.gen_lqr_new(xi_1_node)                
                self.Tree.append(xi_1_node)
                if len(self.Tree) % 10 == 0:
                    print(len(self.Tree))
                self.add_edge(xi, xi_1_node, np.linalg.norm(ui))
                
                if self.near_goal(xi_1):
                    path = self.get_path(xi_1_node)
                    return path

                if len(self.Tree) % 150 == 0:
                    goal_node = Node(self.goal)
                    self.gen_lqr_new(goal_node)
                    nearest_to_goal = self.nearest_neighbor(goal_node)
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
    arm = arm_module.Arm(dt=0.01)

    #state = theta (degrees), theta_dot, time
    
    #Start position 
    start = np.array([0., 0., 0., 0., 0., 0., 0.])

    #end position
    end = np.array([0, 0.5, 0.5, 0, 0.5, 0.5, 3.0]) # This is random...

    #intialize Tree with Xo (starting point)
    rrt = TBRRT(arm, start,end)
    path = rrt.search(10000)