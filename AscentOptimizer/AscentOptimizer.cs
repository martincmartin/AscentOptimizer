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
//		float prevYaw;
//		float prevThrottle;
		FlightCtrlState prevState = new FlightCtrlState();

		SimpleLinearRegression rollToTorque = new SimpleLinearRegression();
		SimpleLinearRegression pitchToTorque = new SimpleLinearRegression();
		SimpleLinearRegression yawToTorque = new SimpleLinearRegression();

		const double regressionHalfLife = 10.0; // In seconds.
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

		public void FixedUpdate()
		{
			Vessel vessel = FlightGlobals.ActiveVessel;
			// TimeWarp.fixedDeltaTime is what *will* be used to compute the *next* velocity, not what
			// was used to compute the current one.  That's in TimeWarp.deltaTime.  Can't they come up with more descriptive
			// names for these things???

			// Vessel.nextVel is always (0, 0, 0), in both FixedUpdate() and Update().
			// Vessel.lastVel is different in every call to FixedUpdate(), as you'd expect with something physics related.

			Vector3d deltaVelocity = FlightGlobals.ActiveVessel.velocityD - prevVelocity;

			print("fixedDT: " + toStr(TimeWarp.fixedDeltaTime) + ", deltaT: " + toStr(TimeWarp.deltaTime) + "\n" +
			      "accel: " + toStr(vessel.acceleration_immediate, 5) + ", from deltaT: " + toStr(deltaVelocity / TimeWarp.deltaTime, 5) + "\n" +
			      "vel: "+toStr(vessel.velocityD) + ", last correction: " + toStr(Krakensbane.GetLastCorrection()) +
			      ", GetFrameVelocity(): "+toStr(Krakensbane.GetFrameVelocity()));
			//////////  Order of calling
			//
			//  So postAutoPilot() and fly() are actually called BEFORE FixedUpdate().
			//////////  Throttle
			//
			// Learning: FlightInputHandler.state.mainThrottle should be used with accel computed in this frame.

			//////////  Yaw
			//
			// Learning: FlightInputHandler.state.yaw should be used with acccel from *future* frames.
			// Or you could use ModuleGimbal, but that's not guaranteed to exist.
			/*
			foreach (var p in vessel.Parts)
			{
//				print("name: " + p.name + ", partname: " + p.partName + ", classname: " + p.ClassName);
//				string mods = "";
				foreach (PartModule m in p.Modules)
				{
					//					mods += "     name: " + m.name + ", GUIName: " + m.GUIName + ", classname: " + m.ClassName + ", modulename: " + m.moduleName + "\n";

					if (m is ModuleGimbal)
					{
						// So actuation (and actuationLocal) change much closer in time to what the physics sees.
						// However, they still change instantaneously and the acceleration still has lag.
						var g = (ModuleGimbal)m;
						// In actuation, z = yaw.
						print("actuation: " + toStr(g.actuation) + ", local: " + toStr(g.actuationLocal) +
							 ", xMult: " + toStr(g.xMult) + ", yMult: " + toStr(g.yMult));
					}
					else if (m is ModuleEngines) {
						var e = (ModuleEngines)m;
						print("requestedThrottle: " + toStr(e.requestedThrottle) + ", currentThrottle: " + toStr(e.currentThrottle) + ", responseRate: " + toStr(e.throttleResponseRate));
					}
				}
//				print(mods);
			}
			*/

			/*
			print("deltaTime: " + TimeWarp.deltaTime + "\n" +
				"vessel.acceleration: " + toStr(vessel.acceleration) +
				  ", from srf_velocity: " + toStr((vessel.srf_velocity - prevVel) / TimeWarp.deltaTime) + "\n"+
				"lastVel: " + toStr(vessel.lastVel) + ", nextVel: " + toStr(vessel.nextVel)+"\n"+
				"orbit.pos: " + toStr(vessel.orbit.pos)+", angularVelocity: " + toStr(vessel.angularVelocityD));
*/


			// Vector3 velWRTVessel = vessel.transform.InverseTransformDirection(vessel.srf_vel_direction);
			// So, what's the formula we need?  A fixed yaw means the thrust is at a fixed angle w.r.t.
			// heading, but it won't (in general) be going through the center of mass, so it will
			// cause some angular acceleration, d^2 theta / d t^2.  So I think we have a second order
			// system after all: d^2 theta / d t^2 = k1 yaw + k2 * d theta / dt.  Although maybe that's it?
			// Maybe there's no term proportional to theta, so its just first order in d theta / dt?
			//
			// It turns out, vessel has both angularVelocity and angularMomentum.  It also has MOI which I assume
			// is Moment of Inertia.

			prevAngularMom = vessel.transform.TransformDirection(vessel.angularMomentum);
			prevVelocity = vessel.velocityD;
			print("------- End FixedUpdate ----------");
		}

		public void Update()
		{
//			Vessel vessel = FlightGlobals.ActiveVessel;

//			print("============= Update ==========\n");
/*			      "lastVel: " + toStr(vessel.lastVel) + ", nextVel: " + toStr(vessel.nextVel));
			      */

			if (GameSettings.MODIFIER_KEY.GetKey() && Input.GetKeyDown(key))
			{
				controlEnabled = !controlEnabled;
				if (controlEnabled)
				{
					FlightGlobals.ActiveVessel.OnFlyByWire += new FlightInputCallback(fly);
					FlightGlobals.ActiveVessel.OnPostAutopilotUpdate += new FlightInputCallback(postAutopilot);
				}
				else {
					FlightGlobals.ActiveVessel.OnFlyByWire -= new FlightInputCallback(fly);
					FlightGlobals.ActiveVessel.OnPostAutopilotUpdate -= new FlightInputCallback(postAutopilot);
				}
			}
		}

		private void postAutopilot(FlightCtrlState s)
		{
			// This seems to affect, not the next FixedUpdate(), but the one after that, or even the one after that!
			// Actually, there seems to be a little bit of hidden state: the yaw control is changing something else
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

			Vessel vessel = FlightGlobals.ActiveVessel;
			Vector3d linearAcc = (vessel.velocityD - prevVelocity) / TimeWarp.deltaTime;

			Vector3 globalAngularMom = vessel.transform.TransformDirection(vessel.angularMomentum);

			// For torque: x: pitch, y: roll, z: yaw
			Vector3 torque = vessel.transform.InverseTransformDirection((globalAngularMom - prevAngularMom) / TimeWarp.deltaTime);
			print("postAutopilot: throttle: " + s.mainThrottle + ", yaw: " + s.yaw+"\n"+
			      "throttle: " + s.mainThrottle + " " + FlightInputHandler.state.mainThrottle+" "+prevState.mainThrottle+" -> " + toStr(linearAcc.magnitude) + ", yaw: " + prevState.yaw+" -> " + toStr(torque.z, 5)+"\n"+
			      "pitch: "+toStr(prevState.pitch)+", roll: "+toStr(prevState.roll)+", torque: "+toStr(torque));//+", transformed: " + toStr(vessel.transform.TransformDirection(vessel.angularMomentum)));

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

			prevState.CopyFrom(s);
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
				FlightGlobals.ActiveVessel.OnFlyByWire -= new FlightInputCallback(fly);
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



			GUILayout.EndVertical();
			GUI.DragWindow();
		}

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
		//
		// 
		//


		private void fly(FlightCtrlState s)
		{
			var vessel = FlightGlobals.ActiveVessel;

			Orbit orbit = FlightGlobals.ActiveVessel.orbit;

			// MOI: y is for roll.
			print("MOI: " + toStr(vessel.MOI));

			// A first goal could just be to kill velocity, i.e. angular momentum.  Could just stay in momentum/torque
			// world, wouldn't need to use MOI to translate into velocity/position.


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
		}
	}
}
