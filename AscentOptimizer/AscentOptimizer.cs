using System;
using UnityEngine;

namespace AscentOptimizer
{
	[KSPAddon(KSPAddon.Startup.Flight, false)]
	public class AscentOptimizer : MonoBehaviour
	{
		public KeyCode key = KeyCode.K;
		public static bool controlEnabled = false;

		public void Start()
		{
			enabled = true;
			// Some day we might add settings, like AeroGUI
			// ConfigNode settings = GameDatabase.Instance.GetConfigNodes("AEROGUI")[0];
		}

		public void Update()
		{
			print("In Update");
			if (GameSettings.MODIFIER_KEY.GetKey() && Input.GetKeyDown(key))
			{
				print("Saw key");
				controlEnabled = !controlEnabled;
				if (controlEnabled)
				{
					print("Turning on");
					FlightGlobals.ActiveVessel.OnFlyByWire += new FlightInputCallback(fly);
				}
				else {
					print("Turning off");
					FlightGlobals.ActiveVessel.OnFlyByWire -= new FlightInputCallback(fly);
				}
			}

		}

		private void fly(FlightCtrlState s)
		{
			print("In fly");
			s.yaw = -0.2f;
			// see http://wiki.kerbalspaceprogram.com/wiki/Module_code_examples

			// see also http://forum.kerbalspaceprogram.com/index.php?/topic/47341-flightctrlstate-not-working-as-i-expect/
		}
	}
}
