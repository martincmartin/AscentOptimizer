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

		SimpleLinearRegression pointing_from_yaw;
		public double prevTime = -1;
		public const double regressionHalfLife = 1.0; // In seconds.
		public double regressionFactor = Math.Log(2) / regressionHalfLife;

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

		public void AddLabel(string text, GUIStyle style = null)
		{
			if (style == null)
			{
				style = defaultLabelStyle;
			}
			GUILayout.BeginHorizontal();
			GUILayout.Label(text, style);
			GUILayout.EndHorizontal();
		}

		public void DrawCloseButton()
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

		public static String toStr(Vector3 v)
		{
			return "(" + v.x.ToString("N3") + ", " + v.y.ToString("N3") + ", " + v.z.ToString("N3") + ")";
		}

		public void DrawWindow(int windowID)
		{
			GUILayout.BeginVertical();
			GUILayout.FlexibleSpace();

			DrawCloseButton();

			Vessel vessel = FlightGlobals.ActiveVessel;

			Orbit orbit = vessel.orbit;
			AddLabel("orbit.an: " + orbit.an.ToString("N3"));
			AddLabel("orbit.altitude: " + orbit.altitude.ToString("N3"));
			// AeroGUI uses nVel = v.srf_velocity.normalized
			// and: angle of attack = Vector3d.Dot(v.transform.forward, nVel);
			// and: side slip = Vector3d.Dot(v.transform.up, Vector3d.Exclude(v.transform.forward, nVel).normalized);
			// and: climb rate = Vector3d.Dot(v.srf_velocity, v.upAxis);
			AddLabel("srf_velocity: " + vessel.srf_velocity);
			AddLabel("srf_vel_direction: " + vessel.srf_vel_direction);
			// "y" is forward.  "x" is yaw (left/right).  "z" is pitch.
			AddLabel("InverseTransformDirection: " + toStr(vessel.transform.InverseTransformDirection(vessel.srf_vel_direction)));

			// vessel.up and vessel.upAxis seem to be the same, a vector pointing up, away
			// from the center of the reference body through the vessel.  Or something.
			AddLabel("vessel.upAxis: " + vessel.upAxis);
			AddLabel("vessel.up: " + vessel.up);

			// Note that orbit.vel and orbit.pos do NOT change as we're sitting on
			// the launch pad with Kerbin rotating.  So the velocity is relative to
			// a rotating frame, but does not subtract out the rotation of that frame.
			AddLabel("orbit.vel: " + orbit.vel);
			AddLabel("orbit.pos: " + orbit.pos);
			AddLabel("getRFrmVel: " + orbit.referenceBody.getRFrmVel(orbit.pos));
			AddLabel("ref bod angular vel: " + orbit.referenceBody.angularVelocity);
			AddLabel("ref bod angular vel Z up: " + orbit.referenceBody.zUpAngularVelocity);

			AddLabel("calculated vel: " + (orbit.vel - orbit.referenceBody.getRFrmVel(orbit.pos)));
			// I think vessel.GetTransform() is the position + orientation of the vessel
			// relative to its start position on the launch pad.  y = "up", out of the plane
			// of planets; +ve = north pole.  I don't understand X vs Z though.
			AddLabel("position from transform: " + vessel.GetTransform().position);
			AddLabel("-----------------------");

			double curTime = FlightGlobals.ActiveVessel.missionTime;
			AddLabel("mission time: " + curTime.ToString("N3"));

			double weight = 1;
			if (prevTime >= 0)
			{
				double prevWeight = Math.Exp(-regressionFactor * (curTime - prevTime));
				pointing_from_yaw.decay(prevWeight);
				weight = 1 - prevWeight;
			}
			pointing_from_yaw.observe(FlightInputHandler.state.yaw, angle, weight);

			prevTime = curTime;

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
