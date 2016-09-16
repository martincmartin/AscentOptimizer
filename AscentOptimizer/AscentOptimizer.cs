using System;
using UnityEngine;

namespace AscentOptimizer
{
	[KSPAddon(KSPAddon.Startup.Flight, false)]
	public class AscentOptimizer : MonoBehaviour
	{
		public KeyCode key = KeyCode.J;
		public static bool controlEnabled = false;

		public void Start()
		{
			enabled = true;
			// Some day we might add settings, like AeroGUI
			// ConfigNode settings = GameDatabase.Instance.GetConfigNodes("AEROGUI")[0];
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
			}

		}

		private void fly(FlightCtrlState s)
		{
			see http://wiki.kerbalspaceprogram.com/wiki/Module_code_examples

			see also http://forum.kerbalspaceprogram.com/index.php?/topic/47341-flightctrlstate-not-working-as-i-expect/
		}
	}
}
