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

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}


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


// Transform from Frenet s,d coordinates to Cartesian x,y with higer accuracy than getXY
vector<double> my_getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int prev_wp = -1;
	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) )){
		prev_wp++;
	}
	int wp2 = (prev_wp+1)%maps_x.size();

	std::vector<double> X, Y;

	// int = ? where ? represent number of waypoints before the prev_wp
    for(int i = 4; i > 0; i--){
    	double wp = (prev_wp-i)%maps_x.size();
    	X.push_back(maps_x[wp]);
    	Y.push_back(maps_y[wp]);    	
    }

    X.push_back(maps_x[prev_wp]);
    Y.push_back(maps_y[prev_wp]);

    // i > ? where ? represent number of waypoints next to prev_xp
    for(int i = 1 ; i <= 4; i++){
    	double wp = (prev_wp+i)%maps_x.size();
    	X.push_back(maps_x[wp]);
    	Y.push_back(maps_y[wp]);
    }

    // transform from global cordinates to local to the previous waypoint
    double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
    for (int i = 0 ; i < X.size() ; i++){
            //traslation transform
            double x_i = X[i] - maps_x[prev_wp];
            double y_i = Y[i] - maps_y[prev_wp];
            //rotation transform
            X[i] =  x_i*cos(heading) + y_i*sin(heading);
            Y[i] = -x_i*sin(heading) + y_i*cos(heading);
	}

	// create a spline
	tk::spline f;

	// set (X,Y) points to the spline
	f.set_points(X, Y);

	// the arc length we want to find is:
	double seg_s = s - maps_s[prev_wp];

	// integrate to calculate the local coordinate x
	double integral = 100;
	double x_1 = seg_s;
	while(integral > seg_s){
		integral = arc_length(f, 0, x_1);
		x_1 -= 0.001;
	}
	x_1 += 0.001;

	//(x_1,y_1) is the point in local coordinates of the s frenet frame
	double y_1 = f(x_1);

	//calculate the angle of the normal vector (to the right side) at point (x_1,y_1)
	double m = f.deriv(1,x_1);
	double theta = atan(-1.0/m);
	if(m < 0) theta *= -1;

	// calculate the (x_2,y_2) point that is equivalen to (s,d) in local frame
	double x_2 = d * cos(theta) + x_1;
	double y_2 = d * sin(theta) + y_1;

	// rotate back to normal after rotating it earlier
    //rotation transform
    double x =  x_2*cos(-heading) + y_2*sin(-heading);
    double y = -x_2*sin(-heading) + y_2*cos(-heading);
    //traslation transform
    x += maps_x[prev_wp];
   	y += maps_y[prev_wp];

	return {x,y};
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

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
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
function of time t that takes the coefficients of a polynomial.
*/
double poly_eval(vector<double> coefficients, double t){
	double total = 0.0;

	for(int i = 0; i < coefficients.size(); i++){
		total += coefficients[i]*pow(t,i);
	}

	return total;
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

  // start in lane 1
  int lane = 1;

  // Have a reference velocity to target
  double ref_vel = 0.0;

  double s_dot_prev = 0;
  double d_dot_prev = 0;
  double first_s = 0;
  double first_d = 0;
  double real_prev_size = 0;


  h.onMessage([&s_dot_prev,&d_dot_prev,&first_s,&first_d,&real_prev_size,&lane,&ref_vel,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
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

          	/*----------------------------------------------------------------------------------------
			BEGIN: Section of code to cheack where other vehicles are.
          	----------------------------------------------------------------------------------------*/
/*
          	if(prev_size > 0){
          		car_s = end_path_s;
          	}

          	bool too_close = false;

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
          				// ref_vel = 29.5; //mph
          				too_close = true;
          				if(lane > 0){
          					lane = 0;
          				}

          			}
          		}
          	}


          	//
          	if(too_close){
          		ref_vel -= 0.224; //0.224 ~= 5 m/s² that is under 10 requirement
          	}
          	else if(ref_vel < 49.5){
          		ref_vel += 0.224;
          	}
*/

          	/*----------------------------------------------------------------------------------------
			END: Section of code to cheack where other vehicles are.
          	----------------------------------------------------------------------------------------*/         


          	/*----------------------------------------------------------------------------------------
			BEGIN: Stay in a lane using splines
          	----------------------------------------------------------------------------------------*/
/*
          	// Create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
          	// Later we will interpolate these waypoints with a spline and fill it in with more points that control speed
          	vector<double> ptsx;
          	vector<double> ptsy;

          	// reference x,y yaw states
          	// either we will reference the satarting point as where the car is or at the previous paths end point
          	double ref_x = car_x;
          	double ref_y = car_y;
          	double ref_yaw = deg2rad(car_yaw);


          	// if previous path is almost empty, use the car as startinf reference
          	if(prev_size < 2){
          		// Use two points that make the path tangent to the car
          		double prev_car_x = car_x - cos(car_yaw);
          		double prev_car_y = car_y - sin(car_yaw);

          		ptsx.push_back(prev_car_x);
          		ptsx.push_back(car_x);
          		
          		ptsy.push_back(prev_car_y);
          		ptsy.push_back(car_y);

          	}else{
          		// Use the previous path's end point as starting reference 

          		//redefine reference state as previous path end point
          		ref_x = previous_path_x[prev_size-1];
          		ref_y = previous_path_y[prev_size-1];

          		double ref_x_prev = previous_path_x[prev_size-2];
          		double ref_y_prev = previous_path_y[prev_size-2];
          		ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

          		// Use the two points that make the path tangent to the previous path's end point
          		ptsx.push_back(ref_x_prev);
          		ptsx.push_back(ref_x);

          		ptsy.push_back(ref_y_prev);
          		ptsy.push_back(ref_y);

          	}

          	// In Frenet add evenly 30m spaced points ahead of the starting reference
          	vector<double> next_wp0 = getXY(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          	vector<double> next_wp1 = getXY(car_s+60, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          	vector<double> next_wp2 = getXY(car_s+90, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

          	//vector<double> foo = my_getXY(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          	//cout << foo[0] << " " << foo[1] << " = "<< next_wp0[0] << " " << next_wp0[1] << endl;

          	ptsx.push_back(next_wp0[0]);
          	ptsx.push_back(next_wp1[0]);
          	ptsx.push_back(next_wp2[0]);

          	ptsy.push_back(next_wp0[1]);
          	ptsy.push_back(next_wp1[1]);
          	ptsy.push_back(next_wp2[1]);

          	for (int i = 0; i < ptsx.size(); i++){
          		// Shift car reference angle to 0 
          		double shift_x = ptsx[i] - ref_x;
          		double shift_y = ptsy[i] - ref_y;

          		ptsx[i] = shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw);
          		ptsy[i] = shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw);

          	}

          	// create a spline
          	tk::spline s;

          	// set (x,y) points to the spline
          	s.set_points(ptsx,ptsy);

          	// Define the actual (x,y) points we will use for the planner
          	// Start with all od the previous path points from last time
          	for(int i = 0; i < previous_path_x.size(); i++){
          		next_x_vals.push_back(previous_path_x[i]);
          		next_y_vals.push_back(previous_path_y[i]);
          	}

          	//Calculate how to break up spline points so that we travel at pur desired reference velocity
          	double target_x = 30.0;
          	double target_y = s(target_x);
          	//double target_dist = sqrt(target_x*target_x + target_y*target_y);
          	double target_dist = sqrt((target_x)*(target_x) + (target_y)*(target_y)); // parentesis

          	double x_add_on = 0; 

          	// Fill up the rest of our path planner after filling it with previous points, here we eill always output 50 points
          	for(int i = 1; i <= 50-previous_path_x.size(); i++){
          		double N = target_dist / (0.02*ref_vel/2.24);
          		double x_point = x_add_on + target_x/N;
          		double y_point = s(x_point);

          		x_add_on = x_point;

          		double x_ref = x_point;
          		double y_ref = y_point;

          		// rotate back to normal after rotating it earlier
          		x_point = x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw);
          		y_point = x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw);

          		x_point += ref_x;
          		y_point += ref_y;

          		next_x_vals.push_back(x_point);
          		next_y_vals.push_back(y_point);

          	}
*/
          	/*----------------------------------------------------------------------------------------
			END: Stay in a lane using splines
          	----------------------------------------------------------------------------------------*/


			/*----------------------------------------------------------------------------------------
			BEGIN: Stay in a lane without splines
          	----------------------------------------------------------------------------------------*/

/*
		    for(int i = 0; i < previous_path_x.size() ; i++){
          		next_x_vals.push_back(previous_path_x[i]);
          		next_y_vals.push_back(previous_path_y[i]);
          	}


          	vector<double> last_frenet;

          	if(prev_size >= 2){
          		int index_1 = prev_size -1;
          		int index_0 = index_1 - 1;
          		double y1 = previous_path_y[index_1];
          		double y0 = previous_path_y[index_0];
          		
          		double x1 = previous_path_x[index_1];
          		double x0 = previous_path_x[index_0];
          		
          		double theta = atan2(y1-y0, x1-x0);

          		last_frenet = getFrenet(previous_path_x[index_1] , previous_path_y[index_1], theta, map_waypoints_x, map_waypoints_y);
          	}else{
          		last_frenet = {car_s, car_d};
          	}

*/

/*          	
          	double dist_inc = 45./111;
          	int N = 50;
		    for(int i = 0; i < N; i++)
		    {
		    	double s = car_s + dist_inc*(i+1);
		    	double d = 6;

		    	vector<double> XY = my_getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
		        
		        next_x_vals.push_back(XY[0]);
		        next_y_vals.push_back(XY[1]);
		        // next_x_vals.push_back(car_x+(dist_inc*i)*cos(deg2rad(car_yaw)));
		        // next_y_vals.push_back(car_y+(dist_inc*i)*sin(deg2rad(car_yaw)));
		    }
*/		    
		    /*----------------------------------------------------------------------------------------
			END: Stay in a lane without splines
          	----------------------------------------------------------------------------------------*/


		    /*----------------------------------------------------------------------------------------
			BEGIN: Calculate intitial condiction s_dot and s_double_dot
          	----------------------------------------------------------------------------------------*/
          	

          	/*----------------------------------------------------------------------------------------
			END: Stay in a lane without splines
          	----------------------------------------------------------------------------------------*/
          	// convert car_speed to s_dot and d_dot

          	int n_prev = real_prev_size - prev_size ;

          	double s_dot = 0;
          	double d_dot = 0;
          	double s_double_dot = 0;
          	double d_double_dot = 0;

          	if( n_prev > 0){
          		s_dot = (car_s - first_s) / 0.02;
          		d_dot = (car_d - first_d) / 0.02;
          	}

          	/*if( n_prev > 1){
          		s_double_dot = (s_dot - s_dot_prev) / 0.04;
          		d_double_dot = (d_dot - d_dot_prev) / 0.04;
          	}

          	s_dot_prev = s_dot;
          	d_dot_prev = d_dot;*/
          	


		    // start location of the vehicle: {s, s_dot, s_double_dot}

          	double T = 1;

          	vector<double> s_start = {car_s, s_dot, s_double_dot};
          	vector <double> s_end = {car_s+10, 15.0/2.24, 0};
          	vector<double> s_trayectory = JMT(s_start, s_end, T);

          	vector<double> d_start = {car_d, d_dot, d_double_dot};
          	vector <double> d_end = {6, 0, 0};
          	vector<double> d_trayectory = JMT(d_start, d_end, T);

          	double t = 0.0;

          	cout <<"Car: "<< car_s << "," << car_d << endl;
          	cout << "Car_speed " << car_speed/2.24 << endl;
          	cout <<"Frenet speeds: "<< s_dot << "," << d_dot <<endl;
          	cout <<"Frenet accc: "<< s_double_dot << "," << d_double_dot <<endl;

          	cout << "***********************\n"; 


      		first_s = car_s;
      		first_d = car_d;

	        t += 0.01;

          	while(t <= T+0.01){
          		double s = poly_eval(s_trayectory, t);
          		double d = poly_eval(d_trayectory, t);

          		cout << s << " , " << d << endl;
          		
          		vector<double> XY = my_getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          		
          		next_x_vals.push_back(XY[0]);
		        next_y_vals.push_back(XY[1]);

		        t += 0.01;
          	}

	        real_prev_size = next_x_vals.size();
		    
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
