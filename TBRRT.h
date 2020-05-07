#include <libfranka.h> // TODO: install
#include <Eigen/Dense> // TODO: install this
#include <vector>
#include <math.h>
#include <random>
#include <limits>

//Lib franka includes
#include <franks/exception.h>
#include <franka/robot.h>

using namespace Eigen;

class TreeNode
{
  private:
    VectorXf state_time;
    TreeNode* parentIdx =  nullptr;
    TreeNode* childIdx = nullptr;
    franka::RobotState& franka_state; //this will have to be saved each time dynamics is applied

  public:
    // --- CONSTRUCTORS --- //
    TreeNode() {}
    TreeNode(VectorXf state_time) : state_time(state_time) {}

    VectorXf get_state() const
    {
      return state_time(seq(0, last-1));
    }

    float get_time() const
    {
      return state_time.tail(1);
    }

    RobotState& getRobotState const
    {
      return franka_state;
    }

};

class TBRRT
{
private:
    std::vector<std::vector<float>> min_state_;
    std::vector<std::vector<float>> max_state_;
    std::vector<TreeNode*> tree;
    std::unordered_map<std::pair<TreeNode*, TreeNode*>, float> edge_costs;
    
    std::array<double, 3> gravity_earth = {0.0, 0.0, -9.8};

    franks::Model model; //should probably be initialized in main and passed into TBRRT when initializing?
    
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

      // min_state = from_franka;
      // max_state = from_franka;

      TreeNode* start_node = new TreeNode(start);
      tree.push_back(start_node);
    }

    // --- MEMBER FUCTIONS --- //
    TreeNode* nearestNeighbor(std::vector<float> xRand)
    {

      /**
       * 1. sample a state, time vector. Maximum time is the time from goal state-time. 
       * 2. consider all treeNodes with a timestamp before the random sample's time as potential nearest neighbors
         3. Compute the norm of distance and velocity between these points. 
         4. Choose the lowest norm as the nearest neighbor. 

            ^^^ end of nearest neighbor

         5. For steering, compute the u between the nearest neighbor and the random neighbor by eqn 12 from the paper. 
            Here, z_doubledot is the finite differenced acceleration, where the dt in the denomiator is the difference in time stamps
            for these respective points

         6. When steering, we will use Franka to apply this U computed above for a time dt (system parameter we choose). 

         7. The point that we end up at after this application of control, will be stored as a new treeNode, with the robotstate from Franka,
            and whose timeStamp is the nearest neighbor's time stamp + dt. 

         8. The point is now added to the tree and we can move on. 

      */



      /**
      * Return the nearest node in the graph to the the random node
      */

      std::vector<float> minState;
      float minControl = std::numeric_limits<double>::infinity();

     //assuming state is a treeNode pointer
      for (auto state : tree)
      {
        franka::RobotState kafka_state = state.getRobotState();

        // https://github.com/frankaemika/libfranka/blob/master/include/franka/robot_state.h
        
        //Get the relevant matrices M, C, and G
        std::array<double, 49> mass_array = model.mass(kafka_state);
        std::array<double, 7> coriolis_array = model.coriolis(kafka_state);
        std::array<double, 7> gravity_array = model.gravity(kafka_state, gravity_earth);
        
        //Convert to Eigen for matrix math
        Eigen::Map<const Eigen::Matrix<double, 7, 7>> M(mass_array.data());
        Eigen::Map<const Eigen::Matrix<double, 7, 1>> C(coriolis_array.data());
        Eigen::Map<const Eigen::Matrix<double, 7, 1>> G(gravity_array.data());
        
        //compute acceleration vector (z_ddot)
        
        //Convert xRand to Eigen vector
        VectirXf rand(xRand.data());
        
        //current state
        VectorXf stateVec = state->get_state();
        
        //accel is finite difference of joint velocities
        VectorXf z_accel = (rand(seq(6,last-1)) - stateVec(seq(6,last-1)))/dt;
        
        //Compute the U from the inertia matrix, coriolois, and gravity matrix from franka
        // u = M * z_accel + C * z_vel + G
        
        VectorXf u =
        
        
        //control = fancy EQN mapping from xRand to state

        if (control < minControl)
        {
            minState = state;
            minControl = control;
        }
      }
    }

    std::vector<float> sample()
    {
      /**
       * Sample the joint space
      */

      std::vector<float> ret;

      // TODO: This might need a generator seed? https://stackoverflow.com/questions/16153589/generating-a-uniform-random-integer-in-c
      for (int i = 0; i < min_state.size(); i++)
      {
          dist = std::uniform_real_distribution(min_state[i], max_state[i]);
          ret.push_back(dist())
      }

      return ret;
    }

    VectorXf steer(VectorXf xi, VectorXf xRand)
    {
      /**
       * Returns ui to move from xi towards xRand
       * :param xi: state representation
       * :param xRand: random state configuration
       */

    }

    VectorXf apply_dynamics(VectorXf xi, VectorXf ui)
    {
      /**
       * :param xi: state_time representation
       * :param ui: 6-dimensional control vector
       */

      // TODO: Return velocities and accelerations for each joint? in state TIME!!!
      // return {theta1_vel, theta2_vel, ...., theta1_acc, theta2_acc... }
      // Are we assuming constant acceleration?

      return; // MAKE SURE to add 1 to the end of the state
    }

    bool is_reachable(std::vector<float> state_time)
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
      float dist_norm = (state_time(seq(0, 5)) - goal(seq(0, 5))).norm();
      float vel_norm = (state_time(seq(6, 11)) - goal(seq(6, 11))).norm();
      float time_norm = (state_time.tail(1) - goal.tail(1)).norm();

      if ((dist_norm < eps[0]) and (vel_norm < eps[1]) and (time_norm < eps[2]))
      {
        return true;
      } else {return false;}
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
        VectorXf ui = steer(xi->get_state(), xRand);
          
        //Note: I think here we should let Franka apply the dynamics and extract the q and dq data for robot_state. Adjusted the line below to reflect this
        
        //should this return an array with robotState and also the state-time vector?
        franka::RobotState xi_1_franka = apply_dynamics(xi->state_time, ui);

        if (is_reachable(xi_1))
        {
          TreeNode* n = new TreeNode(xi_1);
          tree.push_back(n);
          add_edge(xi, xi_1);
          if (in_tolerance(xi_1))
          {
            
          }
        }

      }

      return path;
    }

    // --- DESTRUCTOR --- //
    ~TBRRT()
    {
      for (auto n : tree) {
        delete n;
      }
    }

};
