### Model

The state includes pose of the car (`x`, `y` and orientation `psi`), velocity `v` and error estimates of error in track position `cte` and yaw `epsi`.

The actuators consist of the steering angle `delta` (`-0.43` to `0.43`) and acceleration `a` (`-1` to `1`)

The next state is predicted by updating the state with the incremental state based on current velocity, acceleration, yaw, yaw-error and time-interval `dt` between updates.

### Timestep and Duration

I chose `1.25 s` as the total prediction window and divided that up into timestep `N = 25` and duration `dt = 0.05`. Any wider duration caused wider oscillations in steering, so I settled with this as the final configuration

### Polynomial fitting

Before fitting a polynomial based on the waypoints, they were first translated to compensate for the car's position and rotated by the car's orientation (negative angle to adjust for the simulator). 

This rendered the position of the car to be at origin and at zero orientation.

The polynomial fitted on the transformed waypoints allowed for the `cte` and `epsi` to be calculated.

### Latency

Adding latency to the model simply means that actuations take effect later in the prediction window. So, instead of using previous step's actuation values, we look back `latency_step = latency/dt` steps to find actuation values to use in the update equations. 
