# thermod
<a href="https://github.com/robertladwig/thermod"><img src="images/thermod.png" align="right" height="220" width="220" ></a>

thermod is a simple two-layer water temperature model that assumes that the lake is divided into two volumes: the epilimnion and the hypolimnion. Both layers are divided by a thermocline zone. The entrainment over the thermocline depends on a diffusion coefficient which is a function of the diffusion at neutral stability to the Richardson number.

All equations and derivations are from Steven C. Chapra (2008) 'Surface Water-Quality Modeling' Waveland Press, Inc.

You can run a toy model using the example.R script in `/inst/scripts`.

The package includes example setups for Lough Feeagh (IR) and Lake Mendota (USA).

![](images/2L_visual_result.png)<!-- -->

A comparison between modeled epilimnion (red solid line) and hypolimnion (blue solid line) water temperatures to observed temepratures at 1 (red dots) and 20 m (blue dots) depth, respectively, for Lake Mendota shows that thermod is capable of sufficently replicating real lake dynamics:

![](inst/mendota/2L_compare_mendota.png)<!-- -->

