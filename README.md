# About STRIKE-GOLDD

![GitHub commit activity](https://img.shields.io/github/commit-activity/y/afvillaverde/strike-goldd?label=Commit%20activity&style=plastic)
![GitHub last commit](https://img.shields.io/github/last-commit/afvillaverde/strike-goldd?color=yellow&label=Last%20commit&style=plastic)

STRIKE-GOLDD is a MATLAB toolbox that analyses nonlinear models of ordinary differential equations. It performs a simultaneous assessment of:
- state **observability**,
- parameter structural local **identifiability**,  
- unknown **input observability**. 

Additional features include:
- search for the **Lie symmetries** that underlie non-observability and non-identifiability, 
- automatic **reparameterization** to obtain a fully observable and identifiable model.

Since v4.1.0, it also provides tests for **accessibility** and **controllability**, by integrating the code of the [NLcontrollability](https://github.com/afvillaverde/NLcontrollability/tree/main) toolbox. 

Most of these analyses are performed symbolically, and yield results that are valid for all values of the variables, except for a set of measure zero.

STRIKE-GOLDD was originally created by [Alejandro F. Villaverde](http://afvillaverde.webs.uvigo.gal/), <afvillaverde@uvigo.gal>. It now includes contributions from a number of collaborators.

## Installation and requirements

STRIKE-GOLDD requires a MATLAB installation with the Symbolic Math Toolbox. 

To use STRIKE-GOLDD you just need to:
1. download the code
2. open a MATLAB session
3. run `install.m` from the STRIKE-GOLDD directory
4. define the problem by editing `options.m`
5. run `STRIKE_GOLDD.m`

Alternatively, since v4.0 you can also use it as a Matlab app with a graphical interface:
1. download the code
2. open a MATLAB session
3. right click on `STRIKE_GOLDD_APP.mlapp` and press Run
4. define the problem, method, and options in the pop-up window
5. press the Run button in the pop-up window

If one wishes to use optimization-based decomposition (which is seldom necessary, and currently not recommended), the MATLAB version of the [MEIGO](http://nautilus.iim.csic.es/~gingproc/meigo.html) toolbox is also required.

More information can be found in the [STRIKE-GOLDD manual](STRIKE-GOLDD/doc/STRIKE-GOLDD_manual.pdf).

## Publications

Publication of the methodology (first version of STRIKE-GOLDD):

[Villaverde AF, Barreiro A, Papachristodoulou A (2016). Structural Identifiability of Dynamic Systems Biology Models. *PLoS Computational Biology* 12(10):e1005153, doi:10.1371/journal.pcbi.1005153](http:dx.doi.org/doi:10.1371/journal.pcbi.1005153)

Extension for time-varying inputs (STRIKE-GOLDD 2):

[Villaverde AF, Evans ND, Chappell MJ, Banga JR (2018). Input-dependent structural identifiability of nonlinear systems. *IEEE Control Systems Letters* 3(2):1–6, doi:10.1109/LCSYS.2018.2868608](http://dx.doi.org/doi:10.1109/LCSYS.2018.2868608)

Extension for unknown inputs; FISPO analysis (STRIKE-GOLDD 2.1):

[Villaverde AF, Tsiantis N, Banga JR (2019). Full observability and estimation of unknown inputs, states, and parameters of nonlinear biological models. *Journal of the Royal Society Interface* 16(156), doi:10.1098/rsif.2019.0043](http://dx.doi.org/doi:10.1098/rsif.2019.0043)

Extension for finding Lie symmetries (STRIKE-GOLDD 2.1.6):

[Massonis G & Villaverde AF (2020). Finding and breaking Lie symmetries: implications for structural identifiability and observability in biological modelling. *Symmetry* 12(3), 469, doi:10.3390/sym12030469](https://doi.org/10.3390/sym12030469)

Extension for multi-experiment analysis and implementation of the ORC-DF algorithm (STRIKE-GOLDD 2.2):

[Martínez N & Villaverde AF (2020). Nonlinear observability algorithms with known and unknown inputs: analysis and implementation. *Mathematics* 8(11), 1876, doi:10.3390/math8111876](https://doi.org/10.3390/math8111876)

Extension for automatic reparameterization, AutoRepar (STRIKE-GOLDD 3.0):

[Massonis G, Banga JR, Villaverde AF (2021). AutoRepar: A method to obtain identifiable and observable reparameterizations of dynamic models with mechanistic insight. *International Journal of Robust and Nonlinear Control*, 33(9):5039-5057, doi:10.1002/rnc.5887](https://doi.org/10.1002/rnc.5887)

Extension for graphical interface & new algorithm, ProbObsTest (STRIKE-GOLDD 4.0):

[Díaz-Seoane S, Rey Barreiro X, Villaverde AF (2023). STRIKE-GOLDD 4.0: user-friendly, efficient analysis of structural identifiability and observability. *Bioinformatics*, 39(1), btac748](https://doi.org/10.1093/bioinformatics/btac748)

Description of the accessibility and controllability tests provided with the NLcontrollability tool (STRIKE-GOLDD 4.1):

[Díaz-Seoane S, Barreiro Blas A, Villaverde AF (2023). Controllability and accessibility analysis of nonlinear biosystems. *Computer Methods and Programs in Biomedicine*, 242:107837](https://doi.org/10.1016/j.cmpb.2023.107837)


## Disclaimer

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.
    
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.
