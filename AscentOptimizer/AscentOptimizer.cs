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

		Vector3 prevAngularMom;
		FlightCtrlState prevState = new FlightCtrlState();

		Vector3 deltaAngleUI;
		Vector3 desiredMOMUI;
		Vector3 desiredTorqueUI;

		// Note: a conservative value for prior_x_coeff is a LARGE value.  A large value means we think we only need
		// small control values to produce reasonable torque.
		SimpleLinearRegression rollToTorque = new SimpleLinearRegression(-20.0, 1e-8);
		SimpleLinearRegression pitchToTorque = new SimpleLinearRegression(-20.0, 1e-8);
		SimpleLinearRegression yawToTorque = new SimpleLinearRegression(-20.0, 1e-8);

		const double regressionHalfLife = 5.0; // In seconds.
		double regressionFactor = Math.Log(2) / regressionHalfLife;

		bool doinPhysics = false;
		bool firstPhysicsFrame = true;

		public void Start()
		{
			enabled = true;
			// Some day we might add settings, like AeroGUI
			// ConfigNode settings = GameDatabase.Instance.GetConfigNodes("ASCENTOPT")[0];
			GameEvents.onVesselGoOffRails.Add(VesselGoOffRails);
			GameEvents.onVesselGoOnRails.Add(VesselGoOnRails);

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

		private void VesselGoOffRails(Vessel vessel)
		{
			print("----------  Off Rails  -----------");
			doinPhysics = true;
			firstPhysicsFrame = true;
		}

		private void VesselGoOnRails(Vessel vessel)
		{
			print("----------  On Rails  -----------");
			doinPhysics = false;
		}

		public void Update()
		{
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

		// TimeWarp.fixedDeltaTime is what *will* be used to compute the *next* velocity, not what
		// was used to compute the current one.  That's in TimeWarp.deltaTime.  Can't they come up with more descriptive
		// names for these things???

		// Vessel.nextVel is always (0, 0, 0), in both FixedUpdate() and Update().
		// Vessel.lastVel is different in every call to FixedUpdate(), as you'd expect with something physics related.

		//////////  Throttle
		//
		// FlightInputHandler.state.mainThrottle is the throttle used over the previous time step.

		//////////  Yaw
		//
		// Learning: FlightInputHandler.state.yaw is what will be used over the next time step.


		//////////  Surface vs Orbit
		// vessel.orbit.GetRelativeVel() points prograde relative to orbit.


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

			// Need to guard this against prevState being from a long time ago.  Hmm.
			double prevWeight = Math.Exp(-regressionFactor * TimeWarp.deltaTime);
			double weight = 1 - prevWeight;

			if (Double.IsNaN(prevWeight) || Double.IsNaN(weight))
			{
				print("**********  It's the weight!");
			}

			if (Single.IsNaN(prevState.pitch))
			{
				print("*********  pitch!");
			}
			if (Single.IsNaN(prevState.roll))
			{
				print("*********  roll!");
			}
			if (Single.IsNaN(prevState.yaw))
			{
				print("*********  yaw!");
			}


			if (!Single.IsNaN(torque.x) && !Single.IsNaN(prevState.yaw))
			{
				pitchToTorque.decay(prevWeight);
				pitchToTorque.observe(prevState.pitch, torque.x, weight);
			}

			if (!Single.IsNaN(torque.y) && !Single.IsNaN(prevState.roll))
			{
				rollToTorque.decay(prevWeight);
				rollToTorque.observe(prevState.roll, torque.y, weight);
			}

			if (!Single.IsNaN(torque.z) && !Single.IsNaN(prevState.yaw))
			{
				yawToTorque.decay(prevWeight);
				yawToTorque.observe(prevState.yaw, torque.z, weight);
			}
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

		private void check(float v, String msg)
		{
			if (Single.IsNaN(v))
			{
				print("@@@@@@@@ " + msg);
			}
		}

		private void check(double v, String msg)
		{
			check((float)v, msg);
		}
		
		private void Fly(FlightCtrlState s)
		{
			if (!doinPhysics)
			{
				return;
			}

			if (firstPhysicsFrame)
			{
				print("----------  First frame, skipping update  ----------");
				firstPhysicsFrame = false;
			}
			else {
				UpdateSystemIdentification();
			}

			var vessel = FlightGlobals.ActiveVessel;

			Orbit orbit = vessel.orbit;

			// MOI: y is for roll.
//			print("MOI: " + toStr(vessel.MOI));

			// Note that when torque is zero, angular momentum is constant but angular velcoity can vary.  So,
			// we only have a simple differential equation in yaw/torque/momentum space, not
			// yaw/acceleration/velocity space.  So we multiply everything by moment of inertia and just work there.

			// TODO: Care about the intercept computed by system identification.  For now, we just use the coefficient
			// of the control input.
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
			// Relative to orbit:
			//Vector3d targetDirection = vessel.transform.InverseTransformDirection(vessel.orbit.GetRelativeVel().normalized);
			// Relative to surface
			Vector3d targetDirection = vessel.transform.InverseTransformDirection(vessel.srf_vel_direction);
			check(targetDirection.x, "target direction x");
			check(targetDirection.z, "target direction z");

			Vector3 euler = toEuler(targetDirection);
			check(euler.x, "euler.x");
			check(euler.y, "euler.y");

			// y: yaw
			// x: pitch

			// Translate targetDirection from angle into angular momentum:
			var deltaAngleTimesMOI = new Vector3((float)euler.x * vessel.MOI.x, 0, (float)euler.y * vessel.MOI.z);

			deltaAngleUI = deltaAngleTimesMOI;

			// For angularMomentum: x: pitch, y: roll, z: yaw

			// Let's figure out our desired angular momentum.
			//
			// For now we pretend that all axeses are independent.

			// Since our system identification numbers might be off, and our linear differential equation model is
			// too simple, give ourselves some "breathing room" by planning a trajectory that uses less than the full
			// acceleration we think we can achieve.
			float desiredAccelCoeff = 0.8f;

			// I'm not really sure the best setting for "action time" here.  A lower value means stiffer control.
			// Using TimeWarp.fixedDeltaTime works fine for my simple test rocket, so I use 2 * that just to be safe.
			float actionTime = 2 * TimeWarp.fixedDeltaTime;

			// I originally tried PD control, but since I don't have a desired position trajectory as a function of time,
			// but rather just want to get to the target as quickly as possible by computing a new trajectory every time
			// step, the "p" term was always zero.  And when I tried to add it in anyway, it just caused problems.

			var desiredAngularMomentum =
				new Vector3(getDesiredAngularMomentum(deltaAngleTimesMOI.x, pitchCoeff * desiredAccelCoeff, actionTime),
							0,
				            getDesiredAngularMomentum(deltaAngleTimesMOI.z, yawCoeff * desiredAccelCoeff, actionTime));

			desiredMOMUI = desiredAngularMomentum;


			var desiredTorque = (desiredAngularMomentum - vessel.angularMomentum) / actionTime;
			desiredTorqueUI = desiredTorque;

			check(desiredTorque.x, "desired torque x");
			check(desiredTorque.y, "desired torque y");
			check(desiredTorque.z, "desired torque z");

			check(pitchCoeff, "pitchCoeff");

				
			s.pitch = Math.Max(-1, Math.Min(1, desiredTorque.x / pitchCoeff));
			s.roll = Math.Max(-1, Math.Min(1, desiredTorque.y / rollCoeff));
			s.yaw = Math.Max(-1, Math.Min(1, desiredTorque.z / yawCoeff));


			// "h" is the dot product of position (relative to center of body) and velocity.
			// So, it doesn't take into account the rotation of the body.
			// It's the "specific relative angular momentum."
			//double cosRV = orbit.h.magnitude / orbit.vel.magnitude / orbit.radius;
//			print("an: " + orbit.an);  // One of the anomalies, aka position along orbit?

			// see http://wiki.kerbalspaceprogram.com/wiki/Module_code_examples

			// see also http://forum.kerbalspaceprogram.com/index.php?/topic/47341-flightctrlstate-not-working-as-i-expect/

			prevState.CopyFrom(s);

		}

		public void FixedUpdate()
		{
			Vessel vessel = FlightGlobals.ActiveVessel;

			prevAngularMom = vessel.transform.TransformDirection(vessel.angularMomentum);
		}

	}
}
