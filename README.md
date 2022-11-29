# Psi_StructuralIndicator.jl
Julia script to evaluate the structural indicator Ψ (defined in Foffi &amp; Sciortino, J Phys Chem B 2022).

## How to use
After setting up a working Julia environment, install the required dependencies `DelimitedFiles` and `LightGraphs`
```julia
import Pkg
Pkg.add("DelimitedFiles")
Pkg.add("LightGraphs")
```

Download the sample data files `links.dat` and `distances.dat` and the script `psi.jl` to the same directory.
It will then be sufficient to run (from the Julia REPL)
```julia
include("psi.jl")
```
to reproduce the results in `psi.dat`.

To evaluate Ψ for your own dataset, just process it to match the form of `distances.dat` and `links.dat`,
save it to file with your desired names and then use the function `calc_Ψ`.
The evaluation is extremely simple and can be analyzed from the `psi.jl` script.

## Attribution
If you use this script as a basis for your work including the structural indicator Ψ, please cite our publication

Foffi, R. and Sciortino, F., "Correlated Fluctuations of Structural Indicators Close to The Liquid-Liquid Transition in Supercooled Water", **J. Phys. Chem. B** (...)

## Contact
For any further inquiry, feel free to contact me at rfoffi@ethz.ch
