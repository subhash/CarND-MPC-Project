#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;


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
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

void print_vec(vector<double> v) {
  for (int i=0; i < v.size(); i++) {
    std::cout << v[i] << ", ";
  }
  std::cout << std::endl;
}

void print_vec(Eigen::VectorXd v) {
  for (int i=0; i < v.size(); i++) {
    std::cout << v[i] << ", ";
  }
  std::cout << std::endl;
}

void read_csv(vector<double>& xvals, vector<double>& yvals) {
  std::ifstream file("../lake_track_waypoints.csv");
  std::string header;
  std::getline(file, header);
  while (file) {
    std::string line, xs, ys;
    std::getline(file, line);
    std::stringstream line_stream(line);
    std::getline(line_stream, xs, ',');
    std::getline(line_stream, ys, ',');
    if (!xs.empty() && !ys.empty()) {
      xvals.push_back(std::stod(xs));
      yvals.push_back(std::stod(ys));
    }
  }
  file.close();
}

int main2(){


  vector<double> x_traj, y_traj;
  read_csv(x_traj, y_traj);
  Eigen::VectorXd ptsx = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x_traj.data(), x_traj.size());
  Eigen::VectorXd ptsy = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(y_traj.data(), y_traj.size());
  Eigen::VectorXd coeffs = polyfit(ptsx, ptsy, 3);

//  coeffs << 128.955, 0.338439, -0.00502765;
//  double px = -40.6201, py = 108.73, psi = 3.73367, v = 0.667924;
//  Eigen::VectorXd ptsx(4);
//  ptsx << 30, 40, 50, 60;
//  Eigen::VectorXd ptsy(4);
//  ptsx << -1, -1, -1, -1;
//  coeffs = polyfit(ptsx, ptsy, 1);
//  double px = 100.0, py = -10.0, psi = 0.0, v = 10.0;

  double px = x_traj[0], py = y_traj[0], psi = 0.0, v = 10.0;
  double fx = polyeval(coeffs, px);
  double cte = fx - py;
  double slope = 3*coeffs[3]*px*px + 2*coeffs[2]*px + coeffs[1];
  double psi_ref = atan(slope);
  psi = psi_ref;
  double epsi = psi - psi_ref;
  Eigen::VectorXd state(6);
  state << px, py, psi, v, cte, epsi;


  std::cout << "Coeffs: ";
  print_vec(coeffs);
  std::cout << "Starting state :";
  print_vec(state);
  MPC mpc;

  std::vector<double> x_vals = {state[0]};
  std::vector<double> y_vals = {state[1]};
  std::vector<double> psi_vals = {state[2]};
  std::vector<double> v_vals = {state[3]};
  std::vector<double> cte_vals = {state[4]};
  std::vector<double> epsi_vals = {state[5]};
  std::vector<double> delta_vals = {};
  std::vector<double> a_vals = {};


  double cte_sum = 0;
  int step = 100;
  for (int i=0; i<step; i++) {
    auto vars = mpc.Solve(state, coeffs);
    state << vars[2], vars[3], vars[4], vars[5], vars[6], vars[7];
    std::cout << "State after " << i << " : ";
    print_vec(vars);


    x_vals.push_back(vars[2]);
    y_vals.push_back(vars[3]);
    psi_vals.push_back(vars[4]);
    v_vals.push_back(vars[5]);
    cte_vals.push_back(vars[6]);
    epsi_vals.push_back(vars[7]);

    delta_vals.push_back(vars[0]);
    a_vals.push_back(vars[1]);

    cte_sum += vars[6];

  }

  std::cout << "Average CTE: " << cte_sum/float(step) << std::endl;
  plt::subplot(4, 1, 1);
  plt::title("CTE");
  plt::plot(cte_vals);
  plt::subplot(4, 1, 2);
  plt::title("Delta (Radians)");
  plt::plot(delta_vals);
  plt::subplot(4, 1, 3);
  plt::title("Acceleration");
  plt::plot(a_vals);
  plt::subplot(4, 1, 4);
  plt::title("Velocity");
  plt::plot(v_vals);


  plt::show();

}

int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];

          vector<double> vptsx, vptsy;
          for (int i=0; i<ptsx.size(); i++) {
            double shiftx = ptsx[i] - px, shifty = ptsy[i] - py;
            vptsx.push_back(shiftx * cos(-psi) - shifty * sin(-psi));
            vptsy.push_back(shiftx * sin(-psi) + shifty * cos(-psi));
          }
          double vpx = 0, vpy = 0;
          psi = 0.0;

          Eigen::Vector4d xvals(vptsx.data()), yvals(vptsy.data());
          Eigen::VectorXd coeffs = polyfit(xvals, yvals, 3);
          double fx = coeffs[0];
          double cte = fx - 0;
          double slope = coeffs[1];
          double psi_ref = atan(slope);
          double epsi = 0 - psi_ref;


          Eigen::VectorXd state(6);
          state << vpx, vpy, psi, v, cte, epsi;
          MPC mpc;
          auto results = mpc.Solve(state, coeffs);


          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          double steer_value = -results[0]/deg2rad(25); //-1.0;
          double throttle_value = results[1]; //0.3;

          //std::cout << "CTE: "<< cte << std::endl;


          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory 
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line

          mpc_x_vals = std::vector<double>(results.begin()+8, results.begin()+14);
          mpc_y_vals = std::vector<double>(results.begin()+14, results.begin()+20);

          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          next_x_vals = vptsx;
          next_y_vals = vptsy;

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;


          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          //std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
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
