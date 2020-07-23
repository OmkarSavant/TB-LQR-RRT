# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 19:53:13 2020

@author: samwe
"""
import numpy as np
from rtree import index
import pdb
import networkx as nx
from tqdm import tqdm
import time
import matplotlib.patches as patches
import matplotlib.pyplot as plt

class Obstacle:
    def __init__(self, midpoint, width, height):
        self.mid = midpoint
        self.w = width
        self.h = height
        self.x = midpoint[0] - width/2.0
        self.y = midpoint[1] - height/2.0

# =============================================================================
# class Node:
#     def __init__(self, cost, name):
#         self.name = name
#         self.cost = cost
# 
# class Graph:
#     def __init__(self):
#         self.E = {}
#         
#     def add_edge(child, parent = None):
#         if parent:
#             self.E[child.name] = parent.name
#         else:
#             self.E[child.name] = parent
# =============================================================================

class KDTree():
    def __init__(self, dim=2):
        #TODO: Review this
        
        p = index.Property()
        p.dimension = dim
        self.V = index.Index(interleaved=True, properties=p)  # vertices in an rtree
        self.G = nx.DiGraph()
        
    def add_node(self, val, cost = 0):
        """
        Adds node to KDTree
        """
        val = tuple(val) #TODO: Check this
        self.V.insert(0, val + val, val)
        self.G.add_node(val, cost = cost)
        
    def add_edge(self, child, parent = None):
        self.G.remove_edges_from(self.get_parent(child))
        self.G.add_edge(parent, child)
        
    def get_parent(self, node):
        return self.G.in_edges(node)
        
#%% 

class RRT:
    
    def __init__(self, start, end, obstacles, xlims, ylims):
        """
        start : (start_x, start_y)
        end   : (end_x, end_y)
        obstacles : [ Obstacle_1,  ... ]
        xlims : [x_min, x_max]
        ylims : [y_min, y_max]
        """
        
        # Defined Parameters
        self.R = 0.6 # Robot Radius
        self.G = 0.5 # Goal Radius
        self.max_control = 1
        self.dt = 0.1 # TODO: What should this be?
        self.kRGG = 2*np.e
        
        # Store variables
        self.start = start
        self.end = end
        self.obstacles = obstacles
        self.xlims = xlims
        self.ylims = ylims
        self.min_cost_y = []
        self.min_cost_x = []
        self.best_end = None
        
        # Create KDTree
        self.tree = KDTree()
        
        # Add start node to KDTree
        self.tree.add_node(self.start)
        
        
    def sample(self):
        # returns a state x that is sampled uniformly randomly from the domain
#        return np.random.uniform([self.xlims[0]+self.R, self.ylims[0]+self.R], [self.xlims[1]-self.R, self.ylims[1]-self.R])
        return np.random.uniform([self.xlims[0], self.ylims[0]], [self.xlims[1], self.ylims[1]])

        
        
    def steer(self, x1, x2):
        """
             returns the optimal control trajectory of traveling from x1 to x2
             and the cost
             x1 : nearest
             x2 : rand
             return : x_new, cost
        """
        cost = np.sign(x2 - x1)
        
        return tuple(x2), cost
        
    def is_in_obstacle(self, x):
        # returns True if state x of the robot is incollision with any of the obstacles
        
        for obs in self.obstacles:
            # https://stackoverflow.com/questions/401847/circle-rectangle-collision-detection-intersection
            circleDist_x = abs(x[0] - obs.mid[0])
            circleDist_y = abs(x[1] - obs.mid[1])
            
            if (circleDist_x > (obs.w/2.0 + self.R)):
                continue
            
            if (circleDist_y > (obs.h/2.0 + self.R)):
                continue
            
            if (circleDist_x <= (obs.w/2.0)):
                return True
            
            if (circleDist_y <= (obs.h/2.0)):
                return True
            
            cornerDist_sq = (circleDist_x - obs.w/2.0)**2 + (circleDist_y - obs.h/2.0)**2
            
            if cornerDist_sq <= (self.R**2):
                return True
        
        return False
    
    def is_in_collision(self, start, end):
        """
        Uses interpolation to check for line collision
        """
        
        dist = np.linalg.norm(np.subtract(start, end))
        xs = np.linspace(start, end, int(dist/0.01)) #TODO: Think about parameterized version
        
        for x in xs:
            if self.is_in_obstacle(x):
                return True
        
        return False
        
        
    def nearest(self, x, n=1):
        # finds a node in the tree that is closest to the state x (closest in what metric?)
        x = tuple(x)
        return self.tree.V.nearest(x, num_results=n, objects="raw")
        
#    def connect(self, x):
#        # add state x to the tree
#        # cost of state x = cost of the parent + cost of steer(parent, x)
#        return
#        
#    def rewire(self, x):
#        # rewire all nodes in the tree within the O(gamma (log n/n)ˆ{1/d}} ball
#        # near the state x, update the costs of all rewired neighbors
#        return
    
    def calc_cost(self, x1, x2):
        return np.linalg.norm(np.subtract(x1, x2))
    
    def get_cost(self, node):
        return self.tree.G.nodes[node]["cost"]
    
    def plan(self, n = 10000):
        
        for i in tqdm(range(1, n)):
            x_rand = self.sample()
            x_nearest = next(self.nearest(x_rand))
            
            # x_new, cost = self.steer(x_nearest, x_rand)
            
            x_new = tuple(x_rand)
            
            if not self.is_in_obstacle(x_new):
                X_near = self.nearest(x_new, np.maximum(int(self.kRGG*np.log(i)), 1))
                # self.tree.add_node(x_new)
                x_min = x_nearest
                c_min = self.get_cost(x_nearest) + self.calc_cost(x_nearest, x_rand)
                
                for x_near in X_near:
                    tmp_cost = self.get_cost(x_near) + self.calc_cost(x_near, x_new)
                    if (not self.is_in_collision(x_new, x_near)) and (tmp_cost < c_min):
                        x_min = x_near
                        c_min = tmp_cost
                
                # Can this be simplified
                if not self.is_in_collision(x_min, x_new):
                    self.tree.add_node(x_new, c_min)        
                    self.tree.add_edge(x_new, x_min) # child to parent?
                
                # Rewire
                for x_near in X_near:
                    tmp_cost = self.get_cost(x_new) + self.calc_cost(x_new, x_near)
                    if (not self.is_in_collision(x_new, x_near)) and (tmp_cost < self.get_cost(x_near)):
                        self.tree.add_edge(x_near, x_new)
                        self.tree.G.nodes[x_near]["cost"] = tmp_cost
                        
                if (i % 200 == 0):
                    c_min = np.inf
                    for end in self.nearest(self.end, n=100):
                        if np.linalg.norm(np.subtract(end, self.end)) < (self.R + self.G):
                            tmp_cost = self.get_cost(end)
                            if tmp_cost < c_min:
                                c_min = tmp_cost
                                self.best_end = end
                                
                    if c_min != np.inf:
                        self.min_cost_y.append(c_min)
                        self.min_cost_x.append(i)
                        
                        
                #self.plot()
                
        # Check for path?
        # get nearest to end. Check if intersects
        # self.plot()
        self.optimal_plot()
        self.cost_plot()
                    
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, aspect='equal')
        pos = {key: np.array(key) for key in self.tree.G.nodes()}
        nx.draw_networkx_edges(self.tree.G, pos, ax=ax, width=0.7)
        
        for obs in self.obstacles:
            
            rect = patches.Rectangle((obs.x, obs.y), obs.w, obs.h, facecolor='gray')
            ax.add_patch(rect)
        
        start = patches.Circle(tuple(self.start), 0.2, facecolor='black')
        ax.add_patch(start)
        
        end = patches.Circle(tuple(self.end), self.G, facecolor='green')
        ax.add_patch(end)
        
        fig.savefig('hw4q1.png', dpi=1600)
        
    def optimal_plot(self):
        
        path = nx.shortest_path(self.tree.G, source=tuple(self.start), target=self.best_end)
        
        new_graph = nx.DiGraph()
        
        new_graph.add_node(path[0])
        
        for i in range(1, len(path)):
            new_graph.add_node(path[i])
            new_graph.add_edge(path[i-1], path[i])
        
        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111, aspect='equal')
        
        for obs in self.obstacles:
            
            rect = patches.Rectangle((obs.x, obs.y), obs.w, obs.h, facecolor='gray')
            ax2.add_patch(rect)
        
        start = patches.Circle(tuple(self.start), 0.2, facecolor='black')
        ax2.add_patch(start)
        
        end = patches.Circle(tuple(self.end), self.G, facecolor='green')
        ax2.add_patch(end)
        
        end = patches.Circle(tuple(self.best_end), self.R, facecolor='blue')
        ax2.add_patch(end)
        
        pos = {key: np.array(key) for key in new_graph.nodes()}
        nx.draw_networkx_edges(new_graph, pos, ax=ax2, width=0.7)
        
    def cost_plot(self):
        
        plt.figure()
        plt.scatter(self.min_cost_x, self.min_cost_y)
        plt.title("Minimum Cost vs. Iteration")
        plt.xlabel("Iteration Number")
        plt.ylabel("Minimum Cost")
        

            

#%%
        
obs1 = Obstacle([-3, -4.5], 6, 1)
obs2 = Obstacle([4.5, 2.5], 1, 13)
    
rrt = RRT((0, 0), (8, 8), [obs1, obs2], [-10, 10], [-10, 10])
rrt.plan()