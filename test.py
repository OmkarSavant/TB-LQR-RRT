import numpy as np
from tqdm import tqdm
import control
import pdb
import matplotlib.pyplot as plt


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
        
        self.hasDynamics = False
        self.hasLQR = False

        self.state_time = state_time
        self.M = None
        self.C = None
        self.G = None

        #lqr params
        self.K = None
        self.S = None

        self.parent = np.nan
        self.child = np.nan
    


class TBRRT():
    def __init__(self, start, goal):
        """
        :param start: robot arm starting position
        :param goal: goal position
        """
        
        self.Tree = []
        self.edge_costs = {} # (parent, child) -> cost

        # Store variables
        self.dt = 0.2

        self.goal = goal
        self.min_state_time = np.array([0, -5*np.pi/6, -5*np.pi/6, -20, -20, -20, 0])
        self.max_state_time = np.array([2*np.pi, 5*np.pi/6, 5*np.pi/6, 20, 20, 20, self.goal[-1]])
        self.eps = np.array([2.8, 0.2, 0.03])

        #store coefficients for the M, C, and G matrices 
        self.g = -9.8

        #mass of each link
        self.m1 = 1.87
        self.m2 = 2.2
        self.m3 = 2.2

        #distance from centroid to corresponding axis
        self.l1 = 0.125
        self.l2 = 0.125
        self.l3 = 0.25

        #link length
        self.L1 = 0.25
        self.L2 = 0.25
        self.L3 = 0.5

        #masses of joints
        self.mc1 = 0.3
        self.mc2 = 0.3
        
        #moments of inertia
        self.J1 = 0.0097395
        self.J2 = 0.011458
        self.J3 = 0.045833
        
        #frictional torque coefficients
        self.c1 = 0.7056
        self.c2 = 0.7056
        self.c3 = 0.7056
        
        self.A11 = self.m1*self.l1**2 + self.J1 + (self.m2 + self.m3 + self.mc1 + self.mc2)*self.L1**2
        self.A12 = (self.m2*self.l2 + (self.m3 + self.mc2)*self.L2)*self.L1
        self.A13 = self.m3 * self.l3 * self.L1

        self.A22 = self.m2*self.l2**2 + self.J2 + (self.m3 + self.mc2)*self.L2**2
        self.A23 = self.m3*self.l3*self.L2
        self.A33 = self.m3*self.l3**2 + self.J3
        
        self.B11 = -(self.c1 + self.c2)
        self.B12 = self.c2 + (self.m2*self.l2 + (self.m3 + self.mc2)*self.L2)*self.L1
        self.B13 = self.m3 * self.l3 * self.L1
        
        self.B21 = self.c2 - (self.m2*self.l1 + (self.m3 + self.mc2)*self.L2)*self.L1
        self.B22 = -(self.c2 + self.c3)
        self.B23 = self.c3 + self.m3*self.l3*self.L2

        self.B32 = self.c3 - self.m3*self.l3*self.L2
        self.B33 = -self.c3

        self.C1 = (self.m1*self.l1 + (self.m2 + self.m3 + self.mc1 + self.mc2)*self.L1)*self.g
        self.C2 = (self.m2*self.l2 + (self.m3 + self.mc2)*self.L2)*self.g
        self.C3 = self.m3*self.l3*self.g

        self.Q = np.eye(6) * 100
        self.R = np.eye(2) * 10

        # Add start node to the tree
        start_node = Node(start)
        self.gen_dynamics(start_node)
        self.gen_lqr(start_node)
        self.Tree.append(start_node)

    def sample(self):
        """
        Sample the joint space
        """
        rand_state_time = np.random.uniform(self.min_state_time, self.max_state_time)
        xRand = Node(rand_state_time)

        self.gen_dynamics(xRand)
        self.gen_lqr(xRand)
        
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

            #find the cost to go from the neighbor to xRand
            v_minus_x = np.matrix((node.state_time[0:-1] - xRand.state_time[0:-1]))

            cost = np.matmul(np.matmul(v_minus_x,xRand.S), v_minus_x.T)

            if cost < minCost:
                minNode = node
                minCost = cost


        return minNode #This should be the state-time vector

    def steer(self, xi, xRand):
        """
        Select ui to move from xi toward xrand
        """
        #use the K matrix from xi, and the difference in state should be (xRand-xi)
        x_del = (xi.state_time[0:-1]- xRand.state_time[0:-1]).reshape(6,1)
        # WE DIDNT DO WRAPPING...
                
        ui = -1 * np.matmul(xi.K, x_del)

        print(ui)
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

    def gen_lqr(self, node):

        """ 
        Create a local estimate of the K and S matrices

        Linearize about x, with u = 0 (as per paper)
        """
        node.hasLQR = True

        control_signal = np.zeros([2,1])
        currState = node.state_time[0:6]

        eps = 0.1

        time_step = 0.1
        
        A = np.zeros([len(currState),len(currState)])
        for i in range(len(currState)):
            x = currState.copy()
            x[i] += eps
            x_inc = self.simulate_system(node, x, control_signal, time_step)
            x = currState.copy()
            x[i] -= eps
            x_dec = self.simulate_system(node, x, control_signal, time_step)
            A[:, i] = (x_inc - x_dec) / (2*eps)

        B = np.zeros([len(currState), len(control_signal)])
        for i in range(len(control_signal)):
            u = control_signal.copy()
            u[i] += eps
            x_inc = self.simulate_system(node, currState, u, time_step)
            u = control_signal.copy()
            u[i] -= eps
            x_dec = self.simulate_system(node, currState, u, time_step)
            B[:,i] = (x_inc - x_dec) / (2*eps)
        
        node.K, node.S, E = control.lqr(A, B, self.Q, self.R)
    
    def gen_dynamics(self, node):

        node.hasDynamics = True
        
        #thetas
        t1, t2, t3 = node.state_time[0], node.state_time[1], node.state_time[2]
        
        #theta dots
        td1, td2, td3 = node.state_time[3], node.state_time[4], node.state_time[5]

        node.M = np.matrix([[self.A11, self.A12 * np.cos(t2 - t1), self.A13 * np.cos(t3 - t1)],
                        [self.A12 * np.cos(t2-t1), self.A22, self.A23 * np.cos(t3 - t2)],
                        [self.A13 * np.cos(t3 - t1), self.A23*np.cos(t3 - t2), self.A33]])

        
        node.C = np.matrix([[self.B11, self.B12*np.sin(t2-t1)*td2, self.B13 * np.sin(t3-t1)*td3],
                        [self.B21*np.sin(t2-t1)*td1, self.B22, self.B23*np.sin(t3-t2)*td3],
                        [-self.B13*np.sin(t3-t1)*td1, self.B32*np.sin(t3-t2)*td2, self.B33]])
        
        node.G = np.matrix([[self.C1*np.sin(t1)],[self.C2*np.sin(t2)],[self.C3*np.sin(t3)]])


    def simulate_system(self, node, state, input, time_step = 0.01): 
        """
        Simulate system. Only works on state, NOT state time
        """
        #velocity component
        x_dot = state[3:6].reshape(3,1)
        # D = np.matrix([-input[0]],[input[0]-input[1]],[input[1]])
        D = np.array([-input[0],input[0]-input[1], input[1]]).reshape(-1, 1)
        
        z_ddot = np.matmul(np.linalg.inv(node.M), (np.matmul(node.C, x_dot) + node.G + D))
        updateVec = np.append(x_dot, z_ddot)

        new_state = state + updateVec * time_step        
        return new_state

    def search(self, max_samples):
        """
        Perform TB-RRT Algorithm
        :param max_samples: Number of samples until termination 
        """
        for k in range(max_samples):
            print(k)

            xRand = self.sample()
            xi = self.nearest_neighbor(xRand)
            xi_1, ui = self.steer(xi, xRand)

            if self.is_reachable(xi_1):
                xi_1 = np.append(xi_1, xi.state_time[-1]+self.dt)
                xi_1_node = Node(xi_1)

                #potentially only need gen_dynamics if we are steering
                self.gen_dynamics(xi_1_node)
                self.gen_lqr(xi_1_node)

                self.Tree.append(xi_1_node)
                self.add_edge(xi, xi_1_node, ui)
                if self.near_goal(xi_1):
                    path = self.get_path(xi_1_node)
                    return path
    
    def test_dynamics(self):
        t1s = []
        t2s = []
        t3s = []
        t1ds = []
        t2ds = []
        t3ds = []

        node = Node(np.array([0.1, 0., 0., 0., 0., 0., 0.]))
        self.gen_dynamics(node)
        self.gen_lqr(node)
        
        for t in range(1000):
            new_state = self.simulate_system(node, node.state_time[0:6], [0.0,0.0], time_step = 0.01)

            t1s.append(new_state[0])
            t2s.append(new_state[1])
            t3s.append(new_state[2])
            t1ds.append(new_state[3])
            t2ds.append(new_state[4])
            t3ds.append(new_state[5])

            node = Node(np.append(new_state, t))
            self.gen_dynamics(node)
            self.gen_lqr(node)
        
        a = plt.plot(t1s)
        plt.savefig('t1')

        plt.figure()
        b = plt.plot(t2s)
        plt.figure()
        plt.plot(t3s)

   


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
    start = np.array([0., 0., 0., 0., 0., 0., 0.])

    #end position
    end = np.array([0, 0.5, 0.5, 0, 0.5, 0.5, 3.0]) # This is random...

    #intialize Tree with Xo (starting point)
    rrt = TBRRT(start,end)
    rrt.test_dynamics()
    #path = rrt.search(10000)