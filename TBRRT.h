#pragma once

#include <libfranka.h> // TODO: install
#include <Eigen/Dense> // TODO: install this
#include <vector>
#include <math.h>
#include <random>
#include <limits>

//Lib franka includes
#include <franks/exception.h>
#include <franka/robot.h>
#include <franka/duration.h>
#include <franka/model.h>
#include <franka/rate_limiting.h>
#include <franka/robot.h>


using namespace Eigen;

class TreeNode
{
  private:
    VectorXf state_time;
    TreeNode* parentIdx =  nullptr;
    TreeNode* childIdx = nullptr;
    franka::RobotState* franka_state = nullptr; //this will have to be saved each time dynamics is applied
    

  public:
    // --- CONSTRUCTORS --- //
    TreeNode() {}
    TreeNode(VectorXf state_time) : state_time(state_time) {}
    TreeNode(VectorXf state_time, franka::RobotState* franka_state) : state_time(state_time), franka_state(franka_state) {}

    VectorXf get_state() const
    {
      return state_time(seq(0, last-1));
    }

    float get_time() const
    {
      return state_time.tail(1);
    }

    RobotState& getRobotState() const
    {
      return franka_state;
    }

    // --- DESTRUCTOR --- //
    ~TreeNode()
    {
      delete franka_state;
    }

};

class TBRRT
{
private:
    std::vector<std::vector<float>> min_state;
    std::vector<std::vector<float>> max_state;
    std::vector<TreeNode*> tree;
    std::unordered_map<std::pair<TreeNode*, TreeNode*>, float> edge_costs;
    std::array<double, 3> gravity_earth = {0.0, 0.0, -9.8};

    franka::Robot robot; //should probably be initialized in main and passed into TBRRT when initializing?
    franka::Model model = robot.loadModel();

    VectorXf goal; // This should be state-time vector
    int DoFs;

    std::vector<float> eps {0.5, 0.5, 0.5}; // {dist_tol, speed_tol, time_tol} Chec, this initialization
    float dt;

public:
    // --- CONSTRUCTOR --- //
    TBRRT(VectorXf start, VectorXf goal) : dt(0.02), goal(goal)
    {
      /**
       * :param start: state_time representation of start. start time should be 0
       * :param goal: state_time representation of goal.
       */

      DoFs = (start.size() - 1) / 2;
      //all the franka parameters are at this link:
      // https://frankaemika.github.io/docs/control_parameters.html?highlight=denavit

      //theta in rad, theta dot in rad/sec
      std::vector<float> min_state = {-2.8973,-1.7628,-2.8973,-3.0718,-2.8973,-0.0175,-2.8973,-2.1750,-2.1750,-2.1750,-2.1750,-2.6100,-2.6100,-2.6100,0.0};
      std::vector<float> max_state = {2.8973,1.7628,2.8973,-0.0698,2.8973,3.7525,2.8973,2.1750,2.1750,2.1750,2.1750,2.6100,2.6100,2.6100,goal.tail(1)};

      TreeNode* start_node = new TreeNode(start);
      tree.push_back(start_node);
    }

    // --- MEMBER FUCTIONS --- //
    TreeNode* nearestNeighbor(VectorXf xRand)
    {

      /**
       * Return the nearest node in the graph to the the random node
       * :param xRand: random state-time vector
       * 
       * XXX 1. sample a state, time vector. Maximum time is the time from goal state-time. 
       * XXX 2. consider all treeNodes with a timestamp before the random sample's time as potential nearest neighbors
         XXX 3. Compute the norm of distance and velocity between these points. 
         XXX 4. Choose the lowest norm as the nearest neighbor. 
      */

      TreeNode* minState;
      float minDistance = std::numeric_limits<double>::infinity(); //TODO: Test min control as the metric

     //assuming state is a treeNode pointer
      for (auto node : tree)
      {
        // Filter nodes that happen in the future
        // TODO: Maybe a better heuristic?
        if (node.state_time.tail(1) >= (xRand.tail(1) - dt)) {
          continue;
        }

        float distance = (node.state_time(seq(0, last-1)) - xRand(seq(0, last-1))).norm(); // This is equally weighting position & velocity norm
        
        if (distance < minDistance)
        {
            minState = node;
            minDistance = distance;
        }
      }

      return minState;
    }

    VectorXf sample()
    {
      /**
       * Sample a random state time vector within the allowable space
       * Returns state-time vector
      */

      VectorXf ret;
      std::mt19937 gen(rd()); // https://stackoverflow.com/questions/16153589/generating-a-uniform-random-integer-in-c
      for (int i = 0; i < min_state.size(); i++)
      {
        std::uniform_real_distribution dist(min_state[i], max_state[i]);
        ret(i) = dist(gen);
      }
      return ret;
    }

    TreeNode* steer(TreeNode* xi_node, VectorXf xRand)
    {
      /**
       * Returns ui to move from xi towards xRand
       * :param xi: state-time vector
       * :param xRand: random state-time vector

        5. For steering, compute the u between the nearest neighbor and the random neighbor by eqn 12 from the paper. 
            Here, z_doubledot is the finite differenced acceleration, where the dt in the denomiator is the difference in time stamps
            for these respective points

         6. When steering, we will use Franka to apply this U computed above for a time dt (system parameter we choose). 

         7. The point that we end up at after this application of control, will be stored as a new treeNode, with the robotstate from Franka,
            and whose timeStamp is the nearest neighbor's time stamp + dt. 

         8. The point is now added to the tree and we can move on. 

       */

      VectorXf xi = xi_node.state_time;

      franka::RobotState* franka_state = node.getRobotState();

      // https://github.com/frankaemika/libfranka/blob/master/include/franka/robot_state.h
      
      //Get the relevant matrices M, C, and G
      std::array<double, 49> mass_array = model.mass(*franka_state);
      std::array<double, 7> coriolis_array = model.coriolis(*franka_state);
      std::array<double, 7> gravity_array = model.gravity(*kafka_state, gravity_earth);
      
      //Convert to Eigen for matrix math
      Eigen::Map<const Eigen::Matrix<double, 7, 7>> M(mass_array.data());
      Eigen::Map<const Eigen::Matrix<double, 7, 1>> C_diag(coriolis_array.data());
      Eigen::Matrix<double, 7, 7> C = C_diag.asDiagonal();
      Eigen::Map<const Eigen::Matrix<double, 7, 1>> G(gravity_array.data());
      
      //accel is finite difference of joint velocities
      float timeDiff = xRand.tail(1) - xi.tail(1);
      
      VectorXf z_accel = (xRand(seq(7,last-1)) - xi(seq(7,last-1)))/timeDiff;

      VectorXf z_vel = xi(seq(7,last-1));
      
      //Compute the U from the inertia matrix, coriolois, and gravity matrix from franka
      // u = M * z_accel + C * z_vel + G
      
      // (7x7) * (7x1) + (7 x 7) * (7*1) + (7*1) --> (7x1)
      VectorXf u = M * z_accel + C * z_vel + G; 
  
      //array for the calculated torque for use with rate limiting function
      std::array<double, 7> tau_d_calculated(u.data);

      //apply rate limiter to the contol. May not be necessary if we are only extracting tau_J_d at the end from each state anyways
      std::array<double, 7> tau_d_rate_limited = franka::limitRate(franka::kMaxTorqueRate, tau_d_calculated, franka_state->tau_J_d);

      // Define callback for the joint torque control loop.
      auto torque_control_callback = [&](const franka::RobotState& state, franka::Duration /*period*/) -> franka::Torques {
              
              //return torque command
              return tau_d_rate_limited;
            }

      //start the real-time control loop
      robot.control(torque_control_callback);

      //record the ending state of the robot and save this into a new treeNode in the tree
      franka::RobotState* endingState = new robot.getState();

      //return this new treeNode* so that an edge can be added between xi and the new TreeNode*
      // endingState.q // This is a std::array<double, 7> of the joint positions
      // endingState.dq // This is a std::array<double, 7> of the joint velocities
      VectorXf end_state_time(14);
      end_state_time << ending_state.q.data(), ending_state.dq.data(), timeDiff;

      return new TreeNode(ending_state_time, endingState);
      
    }

    bool is_reachable(VectorXf state_time)
    {
      /**
       * Returns true if the the state does not result in a self collision, otherwise false
       */

      return true;
    }

    void add_edge(TreeNode* parent, TreeNode* child, float cost)
    {
      parent->child = child;
      child->parent = parent;
      edge_costs[{parent, child}] = cost; //TODO: Check this
    }

    bool in_tolerance(TreeNode* n) const
    {
      // We need some concept of dimensionality to parameterize this
      state_time = n->state_time;

      // TODO: Parameterize the indexing based of DoF
      float dist_norm = (state_time(seq(0, 6) - goal(seq(0, 6))).norm();
      float vel_norm = (state_time(seq(7, 13)) - goal(seq(7, 13))).norm();
      float time_norm = (state_time.tail(1) - goal.tail(1)).norm();

      if ((dist_norm < eps[0]) && (vel_norm < eps[1]) && (time_norm < eps[2]))
      {
        return true;
      } else {return false;}
    }

    std::vector<VectorXf> get_path(TreeNode* end_node)
    {
      // TODO: Should we return the state_time or controls? What is the best for visualization?
      // TODO: Should this be VectorXf?
      std::vector<VectorXf> path;
      TreeNode* node = end_node
      while (node->parent != nullptr)
      {
        path.push_back(node->state_time);
        node = node->parent;
      }

      path.push_back(node->state_time);
      return path.reverse();
    }

    std::vector<std::vector> search(int iterations = 10000)
    {
      /**
      * Performs core TB-RRT algorithm
      */
    
      for (int i = 0; i < iterations; i++)
      {
        
        // What is the best way to handle the state_time math within the TreeNode class?

        VectorXf xRand = sample();
        TreeNode* xi = nearestNeighbor(xRand);
        TreeNode* xi_1 = steer(xi, xRand);

        if (is_reachable(xi_1.state_time))
        {
          tree.push_back(xi_1);
          add_edge(xi, xi_1);
          if (in_tolerance(xi_1))
          {
            std::vector<std::vector<float>> path = get_path(xi_1);
            return path;
          }
        }
      }
    }

    // --- DESTRUCTOR --- //
    ~TBRRT()
    {
      for (auto n : tree) {
        delete n;
      }
    }

};
