using System;
using UnityEngine;

namespace AscentOptimizer
{
	[KSPAddon(KSPAddon.Startup.Flight, false)]
	public class AscentOptimizer : MonoBehaviour
	{
		public KeyCode key = KeyCode.K;
		public static bool controlEnabled = false;
		/*
		public void Start()
		{
			// Some day we might add settings, like AeroGUI
			// ConfigNode settings = GameDatabase.Instance.GetConfigNodes("ASCENTOPT")[0];
		}
		*/

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

		private void fly(FlightCtrlState s)
		{
			Orbit orbit = FlightGlobals.ActiveVessel.orbit;

			// "h" is the dot product of position (relative to center of body) and velocity.
			// So, it doesn't take into account the rotation of the body.
			// It's the "specific relative angular momentum."
			double cosRV = orbit.h.magnitude / orbit.vel.magnitude / orbit.radius;
			print("an: " + orbit.an);  // True Anomaly, aka position along orbit, 
			s.yaw = -0.2f;
			// see http://wiki.kerbalspaceprogram.com/wiki/Module_code_examples

			// see also http://forum.kerbalspaceprogram.com/index.php?/topic/47341-flightctrlstate-not-working-as-i-expect/
		}
	}
}
