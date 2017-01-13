using System;
using UnityEngine;

namespace AscentOptimizer
{
	[KSPAddon(KSPAddon.Startup.Flight, false)]
	public class AscentOptimizer : MonoBehaviour
	{
		public KeyCode key = KeyCode.K;
		public static bool controlEnabled = false;
		public static Rect windowPos = new Rect(200, 100, 0, 0);
		public GUIStyle defaultLabelStyle;

		Vector3 prevVelocity;
		Vector3 prevAngularMom;
		FlightCtrlState prevState = new FlightCtrlState();

		Vector3 deltaAngleUI;
		Vector3 desiredMOMUI;
		Vector3 desiredTorqueUI;

		double dUI;
		double pUI;
		double desiredMOMZUI;
		double desiredMOMZOldUI;
		double deltaAngleTimesMOIzUI;

		Vector3 torqueFromMomentumUI;
		Vector3 torqueFromPositionUI;
		double desiredDeltaYawNextTimestepUI;

		// Note: a conservative value for prior_x_coeff is a LARGE value.  A large value means we think we only need
		// small control values to produce reasonable torque.
		SimpleLinearRegression rollToTorque = new SimpleLinearRegression(-20.0, 1e-8);
		SimpleLinearRegression pitchToTorque = new SimpleLinearRegression(-20.0, 1e-8);
		SimpleLinearRegression yawToTorque = new SimpleLinearRegression(-20.0, 1e-8);

		const double regressionHalfLife = 5.0; // In seconds.
		double regressionFactor = Math.Log(2) / regressionHalfLife;

		public void Start()
		{
			enabled = true;
			// Some day we might add settings, like AeroGUI
			// ConfigNode settings = GameDatabase.Instance.GetConfigNodes("ASCENTOPT")[0];
		}

		public void OnGUI()
		{
			if (controlEnabled)
			{
				// Setting this once on object creation wrecks the whole GUI, I don't understand why.  Creating a new
				// one for each frame works great.
				defaultLabelStyle = new GUIStyle(GUI.skin.label);
				defaultLabelStyle.wordWrap = false;

				windowPos = GUILayout.Window("AscentOptimizer".GetHashCode(), windowPos, DrawWindow, "AscentOptimizer");
			}
		}

		public void Update()
		{
			print("------- Update() ----------");
			if (GameSettings.MODIFIER_KEY.GetKey() && Input.GetKeyDown(key))
			{
				controlEnabled = !controlEnabled;
				if (controlEnabled)
				{
					FlightGlobals.ActiveVessel.OnFlyByWire += new FlightInputCallback(Fly);
//					FlightGlobals.ActiveVessel.OnPostAutopilotUpdate += new FlightInputCallback(PostAutopilot);
				}
				else {
					FlightGlobals.ActiveVessel.OnFlyByWire -= new FlightInputCallback(Fly);
//					FlightGlobals.ActiveVessel.OnPostAutopilotUpdate -= new FlightInputCallback(PostAutopilot);
				}
			}
		}


		void AddLabel(string text, GUIStyle style = null)
		{
			if (style == null)
			{
				style = defaultLabelStyle;
			}
			GUILayout.BeginHorizontal();
			GUILayout.Label(text, style);
			GUILayout.EndHorizontal();
		}

		void DrawCloseButton()
		{
			// Enable closing of the window with "x"
			GUIStyle buttonStyle = new GUIStyle(GUI.skin.button);
			buttonStyle.padding = new RectOffset(5, 5, 3, 0);
			buttonStyle.margin = new RectOffset(1, 1, 1, 1);
			buttonStyle.stretchWidth = false;
			buttonStyle.stretchHeight = false;

			GUILayout.BeginHorizontal();
			if (GUILayout.Button("X", buttonStyle))
			{
				controlEnabled = false;
				FlightGlobals.ActiveVessel.OnFlyByWire -= new FlightInputCallback(Fly);
			}
			GUILayout.EndHorizontal();
		}

		static String toStr(Vector3 v, int precision = 3)
		{
			return "(" + toStr(v.x, precision) + ", " + toStr(v.y, precision) + ", " + toStr(v.z, precision) + ") [" +
				toStr(v.magnitude, precision) + "]";
		}

		static String toStr(double x, int precision = 3)
		{
			return x.ToString("N" + precision);
		}

		void DrawLine(SimpleLinearRegression regression, String xLabel)
		{
			double intercept;
			double x_coefficient;
			double x_var;
			regression.solve(out intercept, out x_coefficient, out x_var);

			AddLabel("torque = " + toStr(x_coefficient, 5) + "*" + xLabel + " + " + toStr(intercept, 5) + ", variance: " + toStr(x_var, 5));
		}

		float coeff(SimpleLinearRegression regression)
		{
			double intercept;
			double x_coefficient;
			double x_var;

			regression.solve(out intercept, out x_coefficient, out x_var);

			return (float)x_coefficient;
		}

		public void DrawWindow(int windowID)
		{
			GUILayout.BeginVertical();
			GUILayout.FlexibleSpace();

			DrawCloseButton();

			Vessel vessel = FlightGlobals.ActiveVessel;

			Orbit orbit = vessel.orbit;

			// So, here's what I've learned about frames, etc.
			//
			// altitude: orbit.altitude.
			//
			// srf_velocity: Velocity relative to the surface / atmosphere.  (KSP has no wind, so they're the same.)
			// Zero on the launch pad.  I'm not sure what the coordiante system is, but maybe it doesn't matter.
			//
			// vessel.GetTransform(): The direction is the direction we're facing.  Ignore the position, and the scale
			// is always 1.
			//
			// vessel.transform.InverseTransformDirection(vessel.srf_vel_direction): Our velocity relative to us.
			// "x" is yaw (left/right).  "y" is forward (prograde speed).  "z" is pitch.
			/*


			// Not sure what orbit.an is, but we don't need it.
			// AddLabel("orbit.an: " + orbit.an.ToString("N3"));
			AddLabel("orbit.altitude: " + toStr(orbit.altitude, 0));
			// AeroGUI uses nVel = v.srf_velocity.normalized
			// and: angle of attack = Vector3d.Dot(v.transform.forward, nVel);
			// and: side slip = Vector3d.Dot(v.transform.up, Vector3d.Exclude(v.transform.forward, nVel).normalized);
			// and: climb rate = Vector3d.Dot(v.srf_velocity, v.upAxis);
			AddLabel("srf_velocity: " + toStr(vessel.srf_velocity, 0));
			AddLabel("srf_vel_direction: " + toStr(vessel.srf_vel_direction));
			// The direction we're heading, relative to how we're facing.
			// "y" is forward.  "x" is yaw (left/right).  "z" is pitch.
			AddLabel("InverseTransformDirection: " + toStr(vessel.transform.InverseTransformDirection(vessel.srf_vel_direction)));

			// vessel.up and vessel.upAxis seem to be the same, a vector pointing up, away
			// from the center of the reference body through the vessel.  Or something.
			AddLabel("vessel.upAxis: " + toStr(vessel.upAxis));
			AddLabel("vessel.up: " + toStr(vessel.up));

			// Note that orbit.vel and orbit.pos do NOT change as we're sitting on
			// the launch pad with Kerbin rotating.  So the velocity is relative to
			// a rotating frame, but does not subtract out the rotation of that frame.
			AddLabel("orbit.vel: " + toStr(orbit.vel));
			AddLabel("orbit.pos: " + toStr(orbit.pos, 0));
//			AddLabel("ref bod angular vel: " + toStr(orbit.referenceBody.angularVelocity));
//			AddLabel("ref bod angular vel Z up: " + toStr(orbit.referenceBody.zUpAngularVelocity));

			AddLabel("-----------------------");

			AddLabel("yaw: " + toStr(FlightInputHandler.state.yaw));
			AddLabel("vessel.angularVelocity: " + toStr(vessel.angularVelocity, 6));
			// vessel.angularVelocity seems to be relative to the vessel:
			// x: pitch (W/S). y: roll (Q/E).  z: yaw (A/D).  Excellent!
			// There is one force that depends on absolute angle: gravity.  And when doing a gravity turn,
			// its probably important to take that into account!
			//
			// So the formula we need is something like:
			// d omega / dt = k1 * yaw + k2 * omega + k3 * angle

//			AddLabel("elapsedTime: " + toStr(elapsedTime, 5));
			AddLabel("angularAcc: " + toStr(angularAcc, 5));
*/
			DrawLine(yawToTorque, "yaw");
			DrawLine(pitchToTorque, "pitch");
			DrawLine(rollToTorque, "roll");

			// Vector3 velWRTVessel = vessel.transform.InverseTransformDirection(vessel.srf_vel_direction);
			// So, what's the formula we need?  A fixed yaw means the thrust is at a fixed angle w.r.t.
			// heading, but it won't (in general) be going through the center of mass, so it will
			// cause some angular acceleration, d^2 theta / d t^2.  So I think we have a second order
			// system after all: d^2 theta / d t^2 = k1 yaw + k2 * d theta / dt.  Although maybe that's it?
			// Maybe there's no term proportional to theta, so its just first order in d theta / dt?
			//
			// It turns out, vessel has both angularVelocity and angularMomentum.  It also has MOI which I assume
			// is Moment of Inertia.

			Vector3d targetDirection = vessel.transform.InverseTransformDirection(vessel.orbit.GetRelativeVel().normalized);
			AddLabel("inv transform: " + toStr(targetDirection));
			// x: yaw, y: 1.00, z: pitch

			AddLabel("euler angles: " + toStr(deltaAngleUI));
			// y: yaw
			// x: pitch

			AddLabel("desired MOM: " + toStr(desiredMOMUI));
			AddLabel(" actual MOM: " + toStr(vessel.angularMomentum));

			AddLabel("desiredTorque: " + toStr(desiredTorqueUI));

			AddLabel("pitch: " + toStr(prevState.pitch) + ", roll: " + toStr(prevState.roll) + ", yaw: " + toStr(prevState.yaw));

			AddLabel("--------------------");
			AddLabel("d = "+toStr(dUI));
			AddLabel("p [d * d / 4] = " + toStr(pUI));
			AddLabel("deltaAngleTimesMOI.z = " + toStr(deltaAngleTimesMOIzUI));
			AddLabel("desiredMomentum.z      = " + toStr(desiredMOMZUI));
			AddLabel("desiredMomentum old .z = " + toStr(desiredMOMZOldUI));
			AddLabel("desiredDeltaYawNextTimestepUI = " + toStr(desiredDeltaYawNextTimestepUI));
			AddLabel("torque from momentum: " + toStr(torqueFromMomentumUI.z));
			AddLabel("torque from position: " + toStr(torqueFromPositionUI.z));

			GUILayout.EndVertical();
			GUI.DragWindow();
		}

		/*
		 * Order of operations in Unity:
		 * 
		 * loop:
		 *     postAutoPilot()
		 *     fly()
		 *     FixedUpdate()
		 *     [physics calculations]
		 * 
		 * Update()
		 */

		/********************  Notes on objects and fields  ************************/

		// There seems to be a little bit of hidden state: the yaw control is changing something else
		// (like the angle of the thruster?) that has some inertia and takes some time to get to its final angle.
		//
		// For throttle, ModuleEngine has these fields:
		//
		// float ModuleEngines.currentThrottle
		// The current internal throttle of the engine, which may be different from the current throttle set by the player if useEngineResponseTime is true.

		// float ModuleEngines.engineAccelerationSpeed
		// How quickly the engine spools up when the user-set throttle is higher than currentThrottle.

		// Each frame, if the user throttle is higher than the engine's currentThrottle, currentThrottle is updated according to the formula

		// currentThrottle += (user throttle - currentThrottle) * engineAccelerationSpeed* dt
		// engineAccelerationSpeed has units of inverse seconds.


		/***********************  Notes on physics, math and control theory  ***************************/

		// So what order of control do we have?
		// Our input affects the derivative of our facing, which affects the
		// derivative of our velocity.  So second order, I think.  So the relevant
		// variables are yaw (our input), the difference between our velocity and heading,
		// and the difference between our actual velocity and desired velocity.
		//
		// To a first approximation, the relationship between facing & actual velocity
		// is simple:
		//
		// v(t + dt) = v(t) + acc(t) * dt
		//
		// We can get the total acceleration (which presumably includes drag & whatnot).  Should
		// be of the form a + b * throttle.  Except that thrust vectoring means we change the
		// direction of the thrust, i.e. the direction of b.
		//
		// So lets just go with a simple formulation.
		//
		// I can download math.net numerics, but I'd like to avoid it.  There's a simple
		// closed form solution for 2 independent variables, so I can focus on that.
		//
		// d facing / d t = a + b0 * yaw
		//
		// So, only one variable!
		//
		// To compute that, we need sum(xi), sum(yi), and sum(xi*yi), as well as n.
		//
		// b = sum(wi * (xi - sum(wi*xi)/wT) * (yi - sum(yi - sum(wi*yi)/wT))
		//         ----------------------
		//         sum(wi(xi - sum(wi*xi)/wT)^2)
		//
		// numerator: sum(wi * (xi * yi - xi * y_ - yi * y_ + x_*y_))
		//          = sum(wi * xi * yi) - sum(wi * xi) * y_ - sum(wi * yi) * x_ + wT * x_*y_
		//          = sum(wi * xi * yi) - 2 * sum(wi * xi) * sum(wi * yi) / wT + wT * sum(wi*xi) * sum(wi*yi) / wT^2
		//          = sum(wi * xi * yi) - sum(wi * xi) * sum(wi * yi) / wT
		//
		// denominator = sum(wi * xi * xi - 2 * wi * xi * x_ + wi * x_^2)
		//             = sum(wi * xi * xi) - 2 * sum(wi * xi) * sum(wi * xi) / wT + wT * sum(wi*xi)^2/wT^2
		//             = sum(wi * xi * xi) - sum(wi * xi)^2 / wT
		//
		// So, I need to keep: sum(wi*xi*yi), sum(wi*xi), sum(wi*yi) and sum(wi).
		//
		// What's more, because they're all proportional to wi, its easy to do exponential weighting!

		private void UpdateSystemIdentification()
		{
			Vessel vessel = FlightGlobals.ActiveVessel;
			//Vector3d linearAcc = (vessel.velocityD - prevVelocity) / TimeWarp.deltaTime;

			Vector3 globalAngularMom = vessel.transform.TransformDirection(vessel.angularMomentum);

			// For torque: x: pitch, y: roll, z: yaw
			Vector3 torque = vessel.transform.InverseTransformDirection((globalAngularMom - prevAngularMom) / TimeWarp.deltaTime);
			/*
			print("postAutopilot: throttle: " + s.mainThrottle + ", yaw: " + s.yaw + "\n" +
				  "throttle: " + prevState.mainThrottle + " -> " + toStr(linearAcc.magnitude) + ", yaw: " + prevState.yaw + " -> " + toStr(torque.z, 5) +
				  ", pitch: " + prevState.pitch + " -> " + toStr(torque.x, 5) + ", roll: " + prevState.roll + " -> " + toStr(torque.y, 5));//+", transformed: " + toStr(vessel.transform.TransformDirection(vessel.angularMomentum)));
				  */

			/*			print("deltaT: " + toStr(TimeWarp.deltaTime) + "\n" +
							  "accel: " + toStr(vessel.acceleration_immediate, 5) + ", from deltaT: " + toStr(deltaVelocity / TimeWarp.deltaTime, 5) + "\n" +
							  "vel: " + toStr(vessel.velocityD) + ", last correction: " + toStr(Krakensbane.GetLastCorrection()) +
							  ", GetFrameVelocity(): " + toStr(Krakensbane.GetFrameVelocity()) + "\n" +
							  "throttle: " + FlightInputHandler.state.mainThrottle + ", yaw: " + FlightInputHandler.state.yaw + ", yawAccel: " + toStr(angularAcc.z, 5));

			*/
			// Need to guard this against prevState being from a long time ago.  Hmm.
			double prevWeight = Math.Exp(-regressionFactor * TimeWarp.deltaTime);
			double weight = 1 - prevWeight;

			yawToTorque.decay(prevWeight);
			pitchToTorque.decay(prevWeight);
			rollToTorque.decay(prevWeight);

			yawToTorque.observe(prevState.yaw, torque.z, weight);
			pitchToTorque.observe(prevState.pitch, torque.x, weight);
			rollToTorque.observe(prevState.roll, torque.y, weight);
		}

		private Vector3 toEuler(Vector3d targetDirection)
		{
			// The "alice" vector is "through the looking glass," i.e. mirrored.  X & Z are swapped.
			Vector3d aliceTarget = new Vector3d();
			aliceTarget.x = targetDirection.x;
			aliceTarget.y = targetDirection.z;
			aliceTarget.z = targetDirection.y;

			Vector3 euler = Quaternion.LookRotation(aliceTarget).eulerAngles;
			if (euler.x > 180)
			{
				euler.x -= 360;
			}
			if (euler.y > 180)
			{
				euler.y -= 360;
			}
			euler *= (float)Math.PI / 180;

			return euler;
		}

		private float sgn(float x)
		{
			if (x < 0)
			{
				return -1;
			}
			else {
				return 1;
			}
		}

		/*
		// "E" means "element wise".
		private Vector3 AbsE(Vector3 v)
		{
			return new Vector3(Math.Abs(v.x), Math.Abs(v.y), Math.Abs(v.z));
		}
		*/

		private float getDesiredAngularMomentum(float posTimesMOI, float accel, float timeStep)
		{
			accel = Math.Abs(accel);
			// We want our position as a function of time to be 0.5 * a * t^2, where we're stopped at time t=0 and
			// the current time is negative.  Solve for time:
			float timeToStop = (float)Math.Sqrt(Math.Abs(posTimesMOI) * 2 / accel);
			if (timeToStop <= 2*timeStep)
			{
				// If we slow down at full accel, we'll overshoot.  Instead, compute delta momentum needed to stop
				// in two time steps.  Using one time step can lead to oscillation (if timeStep == TimeWarp.fixedDeltaTime)
				return -posTimesMOI / (2*timeStep);
			}
			return - sgn(posTimesMOI) * accel * (timeToStop - timeStep);
		}

		private void Fly(FlightCtrlState s)
		{
			UpdateSystemIdentification();

			var vessel = FlightGlobals.ActiveVessel;

			Orbit orbit = vessel.orbit;

			// MOI: y is for roll.
//			print("MOI: " + toStr(vessel.MOI));

			// A first goal could just be to kill velocity, i.e. angular momentum.  Could just stay in momentum/torque
			// world, wouldn't need to use MOI to translate into velocity/position.
			//
			// As for what gains we could use: Given that we know the value of deltaTime, we could just choose a value
			// that would cancel the velocity in n timesteps.  If n = 1, we essentially have bang-bang control.  That's
			// pretty optimistic.  :)  We could choose an n like 3 or 5, at least initially.

			// In theory, I should care about the intercept too.  Maybe later.
			//
			// We have that delta_angularMom / dt = torque = control * coeff(XToTorque).
			//
			// Solve for control:
			//
			// control = torque / coeff(XToTorque).
			float rollCoeff = coeff(rollToTorque);
			float pitchCoeff = coeff(pitchToTorque);
			float yawCoeff = coeff(yawToTorque);

			// Let's face an arbitrary direction, for now prograde:
			Vector3d targetDirection = vessel.transform.InverseTransformDirection(vessel.orbit.GetRelativeVel().normalized);

			Vector3 euler = toEuler(targetDirection);

			// y: yaw
			// x: pitch

			// Let's figure out our desired angular momentum.
			//
			// For now we pretend that all axeses are independent.
			//
			// If we multiply angle by MOI, then we have a simple differential equation in terms of the result.

			// Translate targetDirection from angle into angular momentum:
			var deltaAngleTimesMOI = new Vector3((float)euler.x * vessel.MOI.x, 0, (float)euler.y * vessel.MOI.z);

			deltaAngleUI = deltaAngleTimesMOI;

			// For angularMomentum: x: pitch, y: roll, z: yaw

			// The coefficients above tell us the maximum rate of change of angular momentum.  At maximum breaking,
			// we assume a constant deceleration of coeff.  Assuming we stop at time t = 0 (so current time is negative),
			// we have:
			//
			// mom = - t * coeff
			// angleMOI = - 1/2 * t^2 * coeff
			//
			// We want to determine mom from angleMOI, so solve the second one for t and substitute:
			//
			// t = - sqrt(2 * angleMOI / coeff)
			//
			// Substitue in the first:
			//
			// mom = - sgn(angleMOI) * sqrt(2 * angleMOI / coeff) * coeff
			//     = - sgn(angleMOI) * sqrt(2 * angleMOI * coeff)
			//

			float desiredAccelCoeff = 0.8f;

			// I'm not really sure the best setting for "action time" here.  A lower value means stiffer control.
			// Using TimeWarp.fixedDeltaTime works fine for my simple test rocket, so I use 2 * that just to be safe.
			float actionTime = 2 * TimeWarp.fixedDeltaTime;

			var desiredAngularMomentum =
				new Vector3(getDesiredAngularMomentum(deltaAngleTimesMOI.x, pitchCoeff * desiredAccelCoeff, actionTime),
							0,
				            getDesiredAngularMomentum(deltaAngleTimesMOI.z, yawCoeff * desiredAccelCoeff, actionTime));

			desiredMOMUI = desiredAngularMomentum;
			desiredMOMZOldUI = -sgn(deltaAngleTimesMOI.z) * (float)Math.Sqrt(2 * Math.Abs(deltaAngleTimesMOI.z * yawCoeff * desiredAccelCoeff));
			desiredMOMZUI = desiredAngularMomentum.z;

			// Coefficient of angularMomentum that produces a torque required to eliminate in a handfull of timesteps.
			var d = 1 / actionTime;
			// var desiredTorqueTrajectory = d * (desiredAngularMomentum - vessel.angularMomentum);
			// print(toStr(-d) + " * " + toStr(vessel.angularMomentum) + " - " + toStr(p) + " * " + toStr(deltaAngleTimesMOI) + " = " + toStr(desiredTorque));
			// desiredTorqueUI = desiredTorqueTrajectory;

			// The formula above oscillates when near the target direction.  So if we're close, just use standard PID
			// control.  We define "close" as the control inputs needed for PID control are possible (i.e., between -1 and 1),
			// based only on the magnitudes of the position and velocity.

			// Instead, we should probably compute the desired position * MOI & angular momentum at the next timestep,
			// then use PID control to try to achieve.  If we can get to the position within 1 timestep (i.e. if the distance
			// is less than 1/2 a t^2 for our desired max acceleration), we can just use zero for both and we get PID control.


			// For critical damping, we want d = 2 * sqrt(p), i.e.
			// var p = d * d / 4;

			dUI = d;
			// pUI = p;
			deltaAngleTimesMOIzUI = deltaAngleTimesMOI.z;

			//Vector3 deltaAngleTimesMOIChange = new Vector3(vessel.angularMomentum.x, 0, vessel.angularMomentum.z) * TimeWarp.fixedDeltaTime;

			//var desiredTorque = d * (desiredAngularMomentum - vessel.angularMomentum) + p * (desiredAngleTimesMOI - (deltaAngleTimesMOI + deltaAngleTimesMOIChange));
			var desiredTorque = d * (desiredAngularMomentum - vessel.angularMomentum);// + p * (desiredAngleTimesMOI - deltaAngleTimesMOI);
			desiredTorqueUI = desiredTorque;

			// torqueFromMomentumUI = d * (desiredAngularMomentum - vessel.angularMomentum);
//			torqueFromPositionUI = p * (desiredAngleTimesMOI - deltaAngleTimesMOI);
//			desiredDeltaYawNextTimestepUI = desiredAngleTimesMOI.z - deltaAngleTimesMOI.z;


			/*
			var desiredTorquePIDAbs = AbsE(d * vessel.angularMomentum) + AbsE(p * deltaAngleTimesMOI);
			var desiredTorquePID = -d * vessel.angularMomentum - p * deltaAngleTimesMOI;
			Vector3 desiredTorque = new Vector3();
			if (desiredTorquePIDAbs.x / Math.Abs(pitchCoeff) < 1)
			{
				desiredTorque.x = desiredTorquePID.x;
			}
			else {
				desiredTorque.x = desiredTorqueTrajectory.x;
			}

			if (desiredTorquePIDAbs.y / Math.Abs(rollCoeff) < 1)
			{
				desiredTorque.y = desiredTorquePID.y;
			}
			else {
				desiredTorque.y = desiredTorqueTrajectory.y;
			}

//			if (true)
			if (desiredTorquePIDAbs.z / Math.Abs(yawCoeff) < 1)
			{
				desiredTorque.z = desiredTorquePID.z;
				usingPIDYaw = true;
			}
			else {
				desiredTorque.z = desiredTorqueTrajectory.z;
				usingPIDYaw = false;
			}
			*/
				
			s.pitch = Math.Max(-1, Math.Min(1, desiredTorque.x / pitchCoeff));
			s.roll = Math.Max(-1, Math.Min(1, desiredTorque.y / rollCoeff));
			s.yaw = Math.Max(-1, Math.Min(1, desiredTorque.z / yawCoeff));


			// For critical damping, in a system with equation of motion:
			//
			// x.. = -d * x. - p * x
			//
			// We want d = 2*sqrt(p)
			//
			// Also, when far away and making a big change (e.g. from prograde to retrograde), when do we have to start
			// slowing down?  For a constant, negative acceleration that ends at zero velocity:
			//
			// delta_x = 1/2 * a * t^2
			// delta_v = a * t
			//
			// While accelerating toward our goal, v will be increasing and x decreasing, so we neet to start decelerating
			// (at the latest) when the t needed for both of those cross.  Solving for t^2:
			//
			// t^2 = 2 * delta_x / a
			// t^2 = delta_v^2 / a^2
			//
			// 2 * delta_x = delta_v^2 / a
			// a = delta_v^2 / (2 * delta_x)
			//
			// That doesn't fall out of the PID controller, which switches the sign of the acceleration based on v / x.
			//
			// If d^2 = 4*p, then:
			//
			// a = -d * v - d^2/4 * x
			//
			// which equals zero when (and assuming v and x are opposite signs, as above):
			//
			// d * v = d^2/4 * x
			//
			// v = d/4 * x
			//
			// v^2 = d^2 / 16 * x^2
			//
			// So if we set d^2 / 16 = 2 * a * x ( = p / 4), then we'll switch at the time we want.  Although that
			// depends on x at the crossover point, which depends on specifics of the maneuver, e.g. how far away
			// you are from your destination.  So that's a non-starter.
			//
			// I think we'll use PID when we're close to our target, but not when we're far?  Or maybe I'm overthinking
			// this, and we can just use PID all the time.

			// "h" is the dot product of position (relative to center of body) and velocity.
			// So, it doesn't take into account the rotation of the body.
			// It's the "specific relative angular momentum."
			//double cosRV = orbit.h.magnitude / orbit.vel.magnitude / orbit.radius;
//			print("an: " + orbit.an);  // One of the anomalies, aka position along orbit?

			// To get the motion due to the body's rotation, can we take the cross product
			// of the rocket's position (relative to the center of the body) with the body's
			// rotation vector??
//			s.yaw = -0.2f;
			// see http://wiki.kerbalspaceprogram.com/wiki/Module_code_examples

			// see also http://forum.kerbalspaceprogram.com/index.php?/topic/47341-flightctrlstate-not-working-as-i-expect/

			prevState.CopyFrom(s);

		}

		public void FixedUpdate()
		{
			Vessel vessel = FlightGlobals.ActiveVessel;
			// TimeWarp.fixedDeltaTime is what *will* be used to compute the *next* velocity, not what
			// was used to compute the current one.  That's in TimeWarp.deltaTime.  Can't they come up with more descriptive
			// names for these things???

			// Vessel.nextVel is always (0, 0, 0), in both FixedUpdate() and Update().
			// Vessel.lastVel is different in every call to FixedUpdate(), as you'd expect with something physics related.
			/*
						Vector3d deltaVelocity = FlightGlobals.ActiveVessel.velocityD - prevVelocity;

						print("fixedDT: " + toStr(TimeWarp.fixedDeltaTime) + ", deltaT: " + toStr(TimeWarp.deltaTime) + "\n" +
							  "accel: " + toStr(vessel.acceleration_immediate, 5) + ", from deltaT: " + toStr(deltaVelocity / TimeWarp.deltaTime, 5) + "\n" +
							  "vel: "+toStr(vessel.velocityD) + ", last correction: " + toStr(Krakensbane.GetLastCorrection()) +
							  ", GetFrameVelocity(): "+toStr(Krakensbane.GetFrameVelocity()));
							  */
			//////////  Order of calling
			//
			//  So postAutoPilot() and fly() are actually called BEFORE FixedUpdate().
			//////////  Throttle
			//
			// FlightInputHandler.state.mainThrottle is the throttle used over the previous time step.

			//////////  Yaw
			//
			// Learning: FlightInputHandler.state.yaw is what will be used over the next time step.


			prevAngularMom = vessel.transform.TransformDirection(vessel.angularMomentum);

			// Currently only used for printing.
			prevVelocity = vessel.velocityD;
		}

	}
}
