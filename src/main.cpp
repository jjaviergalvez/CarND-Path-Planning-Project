#include "tictoc.h"
#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include <random>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// BEGIN: C O N S T A N T S definition

const double SPEED_FACTOR = 0.45; // To transform from MPH in cartesian coordinate to m/s in Frenet

const int N_SAMPLES = 5;
const vector<double> SIGMA_S = {10.0, 4.0, 2.0}; // s, s_dot, s_double_dot
const vector<double> SIGMA_D = {1.0, 1.0, 1.0};
const double SIGMA_T = 2.0;
const double MAX_JERK = 10.0; // m/s/s/s
const double EXPECTED_JERK_IN_ONE_SEC = 2; // m/s/s
const double MAX_ACCEL = 10.0; // m/s/s
const double EXPECTED_ACC_IN_ONE_SEC = 1; // m/s in frenet frame
const double SPEED_LIMIT = 49.5 * SPEED_FACTOR; //for the moment this speed corepond in a frenet frame
const double VEHICLE_RADIUS = 1.5; // model vehicle as circle to simplify collision detection

// weights of cost functions
const map<string, double> WEIGHTED_COST_FUNCTIONS = {
	{"exceeds_speed_limit_cost", 	1.0},
	//{"negative_speed_cost", 		10.0},
	//{"time_diff_cost",    			1.0},
    //{"s_diff_cost",       			1.0},
    //{"d_diff_cost",       			1.0},
    //{"efficiency_cost",   			1.0},
    {"max_jerk_cost",     			1.0},
    //{"total_jerk_cost",   			1.0},
    {"collision_cost",    			10.0},
    //{"buffer_cost",       			1.0},
    {"max_accel_cost",    			1.0}
    //{"total_accel_cost",  			1.0}
};

// END: C O N S T A N T S definition

//This splines will be used in getXY() function
tk::spline f_0, f_1, f_2, f_3;


// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
// Idea comes from: https://discussions.udacity.com/t/getxy-using-splines-distance-errors/352364/2?u=galvez
vector<double> getXY(double s, double d)
{	
	double x = f_0(s) + d*f_1(s);
	double y = f_2(s) + d*f_3(s);

	return {x,y};
}


/*
	Helper class. Non-ego vehicles move w/ constant acceleration
*/
class Vehicle {
    
public:

	vector<double> _start_state;

	int _id;

	//	Constructor
	Vehicle(vector<double> car){
		_id = car[0];
		double vx = car[3];
		double vy = car[4];
		double s = car[5];
		double d = car[6];
		
		double vel = sqrt(vx*vx + vy*vy);
		_start_state = { s, vel/2, 0,
						 d,     0, 0 };

	}

	// Return the state at time t
	vector<double> state_in(double t){
		auto ib = _start_state.begin();
		auto ie = _start_state.end();

		vector<double> s(ib, ib+3);
		vector<double> d(ie-3, ie);

		vector<double> state = {
	            s[0] + (s[1] * t) + s[2] * t*t / 2.0,
	            s[1] + s[2] * t,
	            s[2],
	            d[0] + (d[1] * t) + d[2] * t*t / 2.0,
	            d[1] + d[2] * t,
	            d[2],
        };

		return state;
	}

};



// Simpson's rule implementation; ingreal a function f(x) over the interval [a,b]
// For more accurate result, the function being integrated need to be relatively smooth over the interval.
// Base on info from https://en.wikipedia.org/wiki/Simpson%27s_rule
template<class Function>
double arc_length(Function& f, double a, double b)
{
	double c = (a+b)/2.0;
	double d = (b-a)/6.0;

	return d * (f.ds(a) + 4.0*f.ds(c) + f.ds(b));
}


/*
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT 
    an array of length 6, each value corresponding to a coefficent in the polynomial 
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

*/

vector<double> JMT(vector< double> start, vector <double> end, double T)
{

    MatrixXd A(3,3);
	VectorXd b(3);

	double si = start[0];
	double si_dot = start[1];
	double si_double_dot = start[2];

	double sf = end[0];
	double sf_dot = end[1];
	double sf_double_dot = end[2];


	A << pow(T,3)  , pow(T,4)   ,  pow(T,5),
		 3*pow(T,2), 4*pow(T,3) ,  5*pow(T,4),
		 6*T 	   , 12*pow(T,2), 20*pow(T,3);

	b << sf - (si + si_dot*T + 0.5 * si_double_dot * pow(T,2) ), 
		 sf_dot - (si_dot + si_double_dot*T), 
		 sf_double_dot - si_double_dot;

	VectorXd x = A.colPivHouseholderQr().solve(b);

    return {si, si_dot, 0.5*si_double_dot, x[0], x[1], x[2]};
    
}

/*
	Function of time t, that takes the coefficients of a polynomial.
*/
double poly_eval(vector<double> coefficients, double t){
	double total = 0.0;
	for(int i = 0; i < coefficients.size(); i++){
		total += coefficients[i]*pow(t,i);
	}

	return total;
}


/*
    Calculates the derivative of a polynomial and returns
    the corresponding coefficients.
*/
vector<double> differentiate(vector<double> coefficients){
	vector<double> new_coefficients;

	for(int i = 1; i < coefficients.size(); i++){
		new_coefficients.push_back( (i*coefficients[i]) );
	}

	return new_coefficients;
}
    
/*
	Takes the coefficients of a polynomial, calculates it derivative of 
	some order, and evaluate it at time t.	
*/
double poly_deriv_eval(vector<double> coefficients, int order, double t){
	
	assert(order>0);
	
	vector<double> new_coefficients = coefficients;
	
	for(int i = 0; i < order; i++){
		new_coefficients = differentiate(new_coefficients);
	}

	double result = poly_eval(new_coefficients, t);

	return result;
}

/*
	This struct is useful in both following cases:

	1. To store variables related to end conditions of a JMT:
		s = {s, s_dot, s_double_dot}
		d = {d, d_dot, d_double_dot}
		t = T

	2. To store variables related to polynomial coefficients trajectory:
		s = {a_0, a_1, a_2, a_3, a_4, a_5} 
		d = {b_0, a_1, b_2, b_3, b_4, b_5} 
		t = T

*/
struct test_case {	
		vector<double> s;
		vector<double> d;
		double t;
};


/*
	Returns a "perturbed" version of the goal.
*/
vector<vector<double>> perturb_goal(vector<double> goal_s, vector<double> goal_d){
	
	vector<vector<double>> perturb_sd;

	double mu, sig;
	default_random_engine gen;

	vector<double> new_s_goal;
	for(int i=1; i < 3; i++){
		mu = goal_s[i];
		sig = SIGMA_S[i];
		normal_distribution<double> gauss(mu,sig);
		new_s_goal.push_back(gauss(gen));
	}

	perturb_sd.push_back(new_s_goal);

	vector<double> new_d_goal;
	for(int i=1; i < 3; i++){
		mu = goal_d[i];
		sig = SIGMA_D[i];
		normal_distribution<double> gauss(mu,sig);
		new_d_goal.push_back(gauss(gen));
	}

	perturb_sd.push_back(new_d_goal);

	return perturb_sd;

}

/*----------------------------------------------------------------------------------------
	BEGIN: Helpers Functions
----------------------------------------------------------------------------------------*/

/*
    A function that returns a value between 0 and 1 for x in the 
    range [0, infinity] and -1 to 1 for x in the range [-infinity, infinity].

    Useful for cost functions.
*/
double logistic(double x){
	return (2.0 / (1 + exp(-x)) - 1.0);
}


/*
	Function that takes as input the differece of velocity and output 
	a smooth acceleration.

*/
double to_acc(double x){
	double x_0 = 8.0, k = 0.5, L = 1.0;
	//double x_0 = 3.0, k = 1.5, L = 1.0;
	//double x_0 = 5.0, k = 1, L = 1.0;
	return ( L / ( 1. + exp(-k*(x - x_0)) ) );
}

/*
    Calculates the closest distance a particular vehicle during a trajectory.
*/
double nearest_approach(test_case traj, Vehicle vehicle){
	
	double t, cur_s, cur_d, targ_s, targ_d, dist;
	vector<double> targ;
	double closest = 99999;

	for(int i = 0; i < 100; i++){
		t = (float)i / 100 * traj.t;
		cur_s = poly_eval(traj.s, t);
		cur_d = poly_eval(traj.d, t);
		targ = vehicle.state_in(t);
		targ_s = targ[0];
		targ_d = targ[3];
		dist = distance(cur_s, cur_d, targ_s, targ_d);
		if(dist < closest){
			closest = dist;
		}
	}

	return closest;
}

/*
    Calculates the closest distance to any vehicle during a trajectory.
*/
double nearest_approach_to_any_vehicle(test_case traj, vector<Vehicle> vehicles){
	double d, closest = 99999.0;

	for(const auto& v: vehicles){
		d = nearest_approach(traj, v);
		if(d < closest){
			closest = d;
		}
	}

	return closest;
}


/*----------------------------------------------------------------------------------------
	BEGIN: Cost Functions
----------------------------------------------------------------------------------------*/

/*
    Penalizes trajectories that span a duration which is longer or 
    shorter than the duration requested.
*/
double time_diff_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	double t = traj.t;
	double cost = logistic(double(abs(t-T)) / T);

	return cost;
}


/*
    Penalizes trajectories whose s coordinate (and derivatives) 
    differ from the goal.
*/
double s_diff_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){

	double cost = 0.0;

	double s = poly_eval(traj.s, traj.t);
	double s_dot = poly_deriv_eval(traj.s, 1, traj.t);
	double s_ddot = poly_deriv_eval(traj.s, 2, traj.t);
	vector<double> S_traj = {s, s_dot, s_ddot};

	vector<double> S_target = target.s;

	for(int i=1; i < 3; i++){
		double diff = abs(S_traj[i] - S_target[i]);
		cost += logistic(diff/SIGMA_S[i]);
	}

	return cost;
}

/*
	Penalizes trajectories whose d coordinate (and derivatives) 
    differ from the goal.
*/
double d_diff_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	double cost = 0.0;

	double d = poly_eval(traj.d, traj.t);
	double d_dot = poly_deriv_eval(traj.d, 1, traj.t);
	double d_ddot = poly_deriv_eval(traj.d, 2, traj.t);
	vector<double> D_traj = {d, d_dot, d_ddot};

	vector<double> D_target = target.d;

	for(int i=1; i < 3; i++){
		double diff = abs(D_traj[i] - D_target[i]);
		cost += logistic(diff/SIGMA_D[i]);
	}

	return cost;
}

/*
    Binary cost function which penalizes collisions.
*/
double collision_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	double nearest = nearest_approach_to_any_vehicle(traj, predictions);

	if(nearest < 2*VEHICLE_RADIUS)
		return 1.0;
	else
		return 0.0;

}

/*
    Penalizes getting close to other vehicles.
*/
double buffer_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	double nearest = nearest_approach_to_any_vehicle(traj, predictions);
    
    return logistic(2*VEHICLE_RADIUS / nearest);
}


/*
    Penalizes getting out of the drivable area at any point of the trajectory
*/
double stays_on_road_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	double t, d;
	double cost = 0.0;

	for(int i = 0; i < 100; i++){
		t = (double)i / 100 * traj.t;
		d = poly_eval(traj.d, t);
		
		if(d < 0.0 || d > 12){
			cost += d;
		}		
	}

	return cost;
}

/*
    Penalizes exceed the speed limit
*/
double exceeds_speed_limit_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	double t, s_dot, d_dot, speed, max_speed = 0.0;
	double cost = 0.0;

	for(int i = 0; i < 100; i++){
		t = (double)i / 100 * traj.t;
		s_dot = poly_deriv_eval(traj.s, 1, t);
		d_dot = poly_deriv_eval(traj.d, 1, t);

		speed = sqrt(s_dot*s_dot + d_dot*d_dot);
		
		if(speed > max_speed)
			max_speed = speed;
	}

	if(max_speed > SPEED_LIMIT)
		return 1.0;
	else
		return 0.0;

	return cost;
}

/*
    Penalizes negative velocity
*/
double negative_speed_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	double t, s_dot, d_dot, speed, max_speed = 0.0;
	double cost = 0.0;
	bool below_zero = false;

	for(int i = 0; i < 100; i++){
		t = (double)i / 100 * traj.t;
		s_dot = poly_deriv_eval(traj.s, 1, t);
		d_dot = poly_deriv_eval(traj.d, 1, t);

		speed = sqrt(s_dot*s_dot + d_dot*d_dot);
		
		if (speed < 0)
			below_zero = true;
	}

	if(below_zero)
		return 1.0;
	else
		return 0.0;

	return cost;
}

/*
    Rewards high average speeds.
*/
double efficiency_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	//TODO
	return 0.0;
}

double total_accel_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	double t, acc, total_acc = 0.0;
	double dt = double(traj.t) / 100.0;

	for(int i = 0; i < 100; i++){
		t = dt * i;
		acc = poly_deriv_eval(traj.s, 2, t);
		total_acc += abs(acc*dt);
	}

	double acc_per_second = total_acc / T;	
    
    return logistic(acc_per_second / EXPECTED_ACC_IN_ONE_SEC );

}

double max_accel_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){
	double acc, t, max_acc = 0;
	double dt = double(traj.t) / 100.0;

	for(int i = 0; i < 100 ; i++){
		t = dt * i;
		acc = poly_deriv_eval(traj.s, 2, t);
		acc = abs(acc);
		if(acc > max_acc)
			max_acc = acc;
	}

	if(max_acc > MAX_ACCEL)
		return 1.0;
	else
		return 0.0;
}

double max_jerk_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){

	double jerk, t, max_jerk = 0;
	double dt = double(traj.t) / 100.0;

	for(int i = 0; i < 100 ; i++){
		t = dt * i;
		jerk = poly_deriv_eval(traj.s, 3, t);
		jerk = abs(jerk);
		if(jerk > max_jerk)
			max_jerk = jerk;
	}

	if(max_jerk > MAX_JERK)
		return 1.0;
	else
		return 0.0;
}

double total_jerk_cost(test_case traj, test_case target, double T, vector<Vehicle> predictions){

	double t, jerk, total_jerk = 0.0;
	double dt = double(traj.t) / 100.0;

	for(int i = 0; i < 100; i++){
		t = dt * i;
		jerk = poly_deriv_eval(traj.s, 3, t);
		total_jerk += abs(jerk*dt);
	}

	double jerk_per_second = total_jerk / T;	
    
    return logistic(jerk_per_second / EXPECTED_JERK_IN_ONE_SEC );
}

// Reflection of cost functions. 
// Idea from https://stackoverflow.com/questions/19473313/how-to-call-a-function-by-its-name-stdstring-in-c
typedef double (*FnPtr)(test_case, test_case, double, vector<Vehicle>);

map<string, FnPtr> cf = {
	{"exceeds_speed_limit_cost",  	exceeds_speed_limit_cost},
	//{"negative_speed_cost", 		negative_speed_cost},
	//{"time_diff_cost",    		  	time_diff_cost},
    //{"s_diff_cost",       		  	s_diff_cost},
    //{"d_diff_cost",					d_diff_cost},
    //{"efficiency_cost",   			efficiency_cost},
    {"max_jerk_cost",     			max_jerk_cost},
    //{"total_jerk_cost",   			total_jerk_cost},
    {"collision_cost",    			collision_cost},
    //{"buffer_cost",      	 		buffer_cost},
    {"max_accel_cost",    			max_accel_cost}
    //{"total_accel_cost",  			total_accel_cost}
};

double calculate_cost(test_case trajectory, test_case target, double goal_t, vector<Vehicle> predictions, bool verbose){
	double cost = 0.0;

	for (const auto& x: WEIGHTED_COST_FUNCTIONS) {
    	auto fname = x.first;
    	double weight = x.second;
    	double new_cost = weight * cf[fname](trajectory, target, goal_t, predictions);
    	cost += new_cost;
    	if(verbose){
    		cout << "cost for '" << fname << "' is \t" << new_cost << endl;
    	}
  	}

	return cost;
}

test_case min_trajectory_cost(vector<test_case> trajectories, test_case target, double T, vector<Vehicle> predictions){
	double cost, min_cost = 99999;
	test_case tr, best;

	for(int i = 0; i < trajectories.size(); i++){
		tr = trajectories[i];
		cost = calculate_cost(tr, target, T, predictions, false);
		if(cost < min_cost){
			min_cost = cost;
			best = tr;
		}
	}

	return best;
}

/*
    Finds the best trajectory according to WEIGHTED_COST_FUNCTIONS (global).

    arguments:
     start_s - [s, s_dot, s_ddot]

     start_d - [d, d_dot, d_ddot]

     target_vehicle - id of leading vehicle (int) which can be used to retrieve
       that vehicle from the "predictions" dictionary. This is the vehicle that 
       we are setting our trajectory relative to.

     delta - a length 6 array indicating the offset we are aiming for between us
       and the target_vehicle. So if at time 5 the target vehicle will be at 
       [100, 10, 0, 0, 0, 0] and delta is [-10, 0, 0, 4, 0, 0], then our goal 
       state for t = 5 will be [90, 10, 0, 4, 0, 0]. This would correspond to a 
       goal of "follow 10 meters behind and 4 meters to the right of target vehicle"

     T - the desired time at which we will be at the goal (relative to now as t=0)

     predictions - dictionary of {v_id : vehicle }. Each vehicle has a method 
       vehicle.state_in(time) which returns a length 6 array giving that vehicle's
       expected [s, s_dot, s_ddot, d, d_dot, d_ddot] state at that time.

    return:
     (best_s, best_d, best_t) where best_s are the 6 coefficients representing s(t)
     best_d gives coefficients for d(t) and best_t gives duration associated w/ 
     this trajectory.
*/
test_case PTG(vector<double> start_s, vector<double> start_d, vector<double> goal_s, vector<double> goal_d, double T, vector<Vehicle> predictions){
	
	// generate alternative goals
	vector<test_case> all_goals;
	test_case goals;
	double timestep = 0.5;
	double t = T - 4 * timestep;
	
	while(t <= T + 4 * timestep){
		goals.s = goal_s;
		goals.d = goal_d;
		goals.t  = t;

		all_goals.push_back(goals);

		for(int i = 0; i < N_SAMPLES; i++){
			vector<vector<double>> perturbed = perturb_goal(goal_s, goal_d);
			goals.s = perturbed[0];
			goals.d = perturbed[1];
			all_goals.push_back(goals);
		}
		t += timestep;
	}	

	// find best trajectory
	vector<test_case> trajectories;
	test_case trajectory;
	for(int i = 1; i < all_goals.size(); i++){
		vector<double> s_goal = all_goals[i].s;
		vector<double> d_goal = all_goals[i].d;
		double t = all_goals[i].t;

		vector<double> s_coefficients = JMT(start_s, s_goal, t);
		vector<double> d_coefficients = JMT(start_d, d_goal, t);
		
		trajectory.s = s_coefficients;
		trajectory.d = d_coefficients;
		trajectory.t = t;

		trajectories.push_back(trajectory);
	}

	test_case target, best;
	target.s = goal_s;
	target.d = goal_d;
	target.t = T;
	best = min_trajectory_cost(trajectories, target, T, predictions);

	return best;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // set (X,Y) points to the splines that will be used in getXY() function
  f_0.set_points(map_waypoints_s, map_waypoints_x);
  f_1.set_points(map_waypoints_s, map_waypoints_dx);
  f_2.set_points(map_waypoints_s, map_waypoints_y);
  f_3.set_points(map_waypoints_s, map_waypoints_dy);

  // start in lane 1
  double lane = 1;

  // Have a reference velocity to target
  double ref_vel = 47.0 * SPEED_FACTOR;

  vector<double> prev_s_coeff, prev_d_coeff;
  double last_s, last_d;
  double t = 0;

  h.onMessage([&last_s,&last_d,&t,&prev_s_coeff,&prev_d_coeff,&lane,&ref_vel,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	vector<double> previous_path_x = j[1]["previous_path_x"];
          	vector<double> previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

          	int prev_size = previous_path_x.size();

          	if(prev_size > 0){
          		car_s = last_s;
          	}

          	// Looking for if is another car in the same lane in front of us 

          	bool too_close = false;
          	int id_to_follow;
          	// Find ref_v to use
          	for(int i = 0; i < sensor_fusion.size(); i++){
          		// Car is in my lane
          		float d = sensor_fusion[i][6];
          		if(d < (2+4*lane+2) && d > (2+4*lane-2)){
          			double vx = sensor_fusion[i][3];
          			double vy = sensor_fusion[i][4];
          			double check_speed = sqrt(vx*vx + vy*vy);
          			double check_car_s = sensor_fusion[i][5];

          			check_car_s += ((double)prev_size * 0.02 * check_speed); // if using previous points can project s value out
          			// check s values greater than mine and s group
          			if((check_car_s > car_s) && (check_car_s-car_s) < 30){
          				// Do some logic here, lower reference velocity so we dont crach into the car infront of us, could
          				// also flag to try to change lanes.
          				//ref_vel = 29.5; //mph
          				too_close = true;
          				id_to_follow = i;
          			}
          		}
          	}          	

          	// Make predictions base on sensor fusion data
          	vector<Vehicle> predictions;
          	for(int i = 0; i < sensor_fusion.size(); i++){
          		vector<double> car = sensor_fusion[i];
          		Vehicle v(car);
          		predictions.push_back(v);
          	}

          	// Estimate the initial state
          	double s = car_s;
          	double d = car_d;
          	double s_dot = 0;
          	double d_dot = 0;
          	double s_ddot = 0;
          	double d_ddot = 0;
          	if(prev_s_coeff.size() != 0){
          		s = last_s;
          		d = last_d;
          		s_dot = poly_deriv_eval(prev_s_coeff, 1, t);
          		d_dot = poly_deriv_eval(prev_d_coeff, 1, t);
          		s_ddot = poly_deriv_eval(prev_s_coeff, 2, t);
          		d_ddot = poly_deriv_eval(prev_d_coeff, 2, t);
          	}
          	vector<double> s_start = {s, s_dot, s_ddot};
	        vector<double> d_start = {d, d_dot, d_ddot};

	        // Define the end-state
	        vector <double> s_end, d_end;
	        double dist, T, acc, diff_vel;
	        test_case trajectory_to_execute;

	        if(too_close){
	        	cout << "CAR IN FRONT";

	        	double vx = sensor_fusion[id_to_follow][3];
          		double vy = sensor_fusion[id_to_follow][4];
          		double check_speed = sqrt(vx*vx + vy*vy);
          		double check_car_s = sensor_fusion[id_to_follow][5];
          		check_car_s += ((double)prev_size * 0.02 * check_speed);
          		
          		ref_vel = check_speed;
          		dist = check_car_s - car_s;

          		T = 2 * dist / (s_dot + check_speed); //formula 2
          		
        		//cout << "Consider Change Lane" << endl;

        		vector<string> states = {"LCL", "LCR"};
			    if(lane == 0)
			        states.erase(states.begin()); // remove LCL
			    if(lane == 2)
			        states.erase(states.end()); // remove LCR
		    
			    map<string, test_case> state_tr;
			    double min_cost = 99999.;
			    string state_min_cost;
			    double test_lane;

			    for(const auto& state:states){
			    	s_end = {check_car_s, ref_vel, 0};

			    	if(state == "LCL")
			    		test_lane = lane - 1;

			    	if(state == "LCR")
			    		test_lane = lane + 1;

		    		d_end = {2+4*test_lane, 0, 0};
		    		//test_case tr = PTG(s_start, d_start, s_end, d_end, T, predictions);
		    		test_case tr;
		    		tr.s = JMT(s_start, s_end, T);
	        		tr.d = JMT(d_start, d_end, T);
	        		tr.t = T;


		    		test_case target;
					target.s = s_end;
					target.d = d_end;
					target.t = T;

					//cout << state << endl;
					double cost = calculate_cost(tr, target, T, predictions, false);

					if(cost < min_cost){
						min_cost = cost;
						state_min_cost = state;
					}

					state_tr.insert({state, tr});
			    }

			    //if(min_cost < 5){
			    if(false){
			    	//cout << state_min_cost << endl;

			    	if(state_min_cost == "LCL"){
			    		cout << " -> lane change left" << endl;
			    		lane -= 1;
			    	}

			    	if(state_min_cost == "LCR"){
			    		lane += 1;
			    		cout << " -> lane change right" << endl;
			    	}

			    	trajectory_to_execute.s = state_tr[state_min_cost].s;
        			trajectory_to_execute.d = state_tr[state_min_cost].d;
			    }
			    else{
			    	cout << " -> follow the car in front" << endl;
			    	s_end = {check_car_s - 0.5, ref_vel, 0};
			    	d_end = {2+4*lane, 0, 0};
          			trajectory_to_execute.s = JMT(s_start, s_end, T);
	        		trajectory_to_execute.d = JMT(d_start, d_end, T);
			    }

	        	
	        }
	        else{
	        	cout << "FREE";

	        	ref_vel = 47 * SPEED_FACTOR;
	          	diff_vel = abs(ref_vel - s_dot); //use abs to consider little bumpings 
	          	
	          	if(diff_vel < 1){
	          		// cruise control
	          		cout << " -> cruise control" << endl;
	          		T = 1;
	          		dist = ref_vel*T; //formula 3 with cero acc
	          	}
	          	else{ 
	          		// accelerate
	          		cout << "-> speeding up" << endl;
		          	acc = 6 * to_acc(diff_vel);
	          		T = diff_vel/acc; //Formula 1
	          		dist = s_dot*T + T*T*acc/2; //formula 3
          		}	          	
	          	
	          	s_end = {s+dist, ref_vel, 0};
	          	d_end = {2+4*lane, 0, 0};


	          	trajectory_to_execute.s = JMT(s_start, s_end, T);
	        	trajectory_to_execute.d = JMT(d_start, d_end, T);

	        }

	        // Trayectory planner
	        //trajectory_to_execute = PTG(s_start, d_start, s_end, d_end, T, predictions);
	        

          	// Send Values to the controller.
	        // First the path that has not yet follower;
          	next_x_vals = previous_path_x;
	        next_y_vals = previous_path_y;

	        // And then, the values of this trajectory
          	t = 0.0;
          	for(int i = 0; i < 70-previous_path_x.size(); i++){
          		t += 0.02;
          		s = poly_eval(trajectory_to_execute.s, t);
          		d = poly_eval(trajectory_to_execute.d, t);
          		//cout << s << " , " << d << endl;
          		vector<double> XY = getXY(s, d);

          		next_x_vals.push_back(XY[0]);
          		next_y_vals.push_back(XY[1]);
          	}

          	// save some useful variables for the next iteration
          	last_s = s;
          	last_d = d;
          	prev_s_coeff = trajectory_to_execute.s;
	        prev_d_coeff = trajectory_to_execute.d;


		    // TODO END
		    
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
