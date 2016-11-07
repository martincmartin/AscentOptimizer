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

		double prevTime = -10;
		Vector3d prevAngularVel;

		double elapsedTime;
		Vector3d angularAcc;

		SimpleLinearRegression pointing_from_yaw = new SimpleLinearRegression();
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

		public void Update()
		{
			Vessel vessel = FlightGlobals.ActiveVessel;
			double curTime = vessel.missionTime;
			elapsedTime = curTime - prevTime;
			// print("curTime: " + curTime + ", prevTime: " + prevTime + ", elapsed: " + elapsedTime);
			if (elapsedTime <= 1.0 && elapsedTime > 0.0)
			{
				double prevWeight = Math.Exp(-regressionFactor * elapsedTime);
				pointing_from_yaw.decay(prevWeight);
				double weight = 1 - prevWeight;

				angularAcc = (vessel.angularVelocityD - prevAngularVel) / elapsedTime;

				pointing_from_yaw.observe(FlightInputHandler.state.yaw, angularAcc.z, weight);
			}
			// Vector3 velWRTVessel = vessel.transform.InverseTransformDirection(vessel.srf_vel_direction);
			// So, what's the formula we need?  A fixed yaw means the thrust is at a fixed angle w.r.t.
			// heading, but it won't (in general) be going through the center of mass, so it will
			// cause some angular acceleration, d^2 theta / d t^2.  So I think we have a second order
			// system after all: d^2 theta / d t^2 = k1 yaw + k2 * d theta / dt.  Although maybe that's it?
			// Maybe there's no term proportional to theta, so its just first order in d theta / dt?
			//
			// It turns out, vessel has both angularVelocity and angularMomentum.  It also has MOI which I assume
			// is Moment of Inertia.

			prevTime = curTime;
			prevAngularVel = vessel.angularVelocityD;


			if (GameSettings.MODIFIER_KEY.GetKey() && Input.GetKeyDown(key))
			{
				controlEnabled = !controlEnabled;
				if (controlEnabled)
				{
					FlightGlobals.ActiveVessel.OnFlyByWire += new FlightInputCallback(fly);
				}
				else {
					FlightGlobals.ActiveVessel.OnFlyByWire -= new FlightInputCallback(fly);
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
				FlightGlobals.ActiveVessel.OnFlyByWire -= new FlightInputCallback(fly);
			}
			GUILayout.EndHorizontal();
		}

		static String toStr(Vector3 v, int precision = 3)
		{
			return "(" + toStr(v.x, precision) + ", " + toStr(v.y, precision) + ", " + toStr(v.z, precision) + ")";
		}

		static String toStr(double x, int precision = 3)
		{
			return x.ToString("N" + precision);
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

			AddLabel("mission time: " + toStr(vessel.missionTime));
			AddLabel("yaw: " + toStr(FlightInputHandler.state.yaw));
			AddLabel("vessel.angularVelocity: " + toStr(vessel.angularVelocity, 6));
			// vessel.angularVelocity seems to be relative to the vessel:
			// x: pitch (W/S). y: roll (Q/E).  z: yaw (A/D).  Excellent!
			// There is one force that depends on absolute angle: gravity.  And when doing a gravity turn,
			// its probably important to take that into account!
			//
			// So the formula we need is something like:
			// d omega / dt = k1 * yaw + k2 * omega + k3 * angle

			AddLabel("elapsedTime: " + toStr(elapsedTime, 5));
			AddLabel("angularAcc: " + toStr(angularAcc, 5));

			double intercept;
			double x_coefficient;
			double x_var;
			pointing_from_yaw.solve(out intercept, out x_coefficient, out x_var);

			AddLabel("angular accel = " + toStr(x_coefficient, 5) + "*yaw + " + toStr(intercept, 5) + ", x_var: " + toStr(x_var, 5));
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
			Orbit orbit = FlightGlobals.ActiveVessel.orbit;

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
