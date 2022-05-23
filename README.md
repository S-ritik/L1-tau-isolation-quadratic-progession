# L1-tau-isolation-quadratic-progession
Based on https://github.com/sandeepbhowmik1/L1TauIsolation

Workflow the repository

Use MakeIsoLUT codes to build Relaxation and then make WPs

MakeIsoLUT/MakeTau_Iso_LUT.C is used to make lut (code not updated for quadratic proggression)

Then apply calibration ApplyCalibrationZeroBias.C

Use rate_Calculation.C to calaculate rate as the function of pt threshold

Macro/Compare_Rate_Plots.C gives threshold at fix rate

Then use ApplyIsolationForTurnOns.C to out isolation cut using defined options and get turn-ons

Plots turnons using Compare_TurnOn_Plots.C

All the root files created and required are kept in inputfiles folder
