## Sioux Falls Example
<img src="http://www.bgu.ac.il/~bargera/tntp/SiouxFalls/SiouxFallsMap_AAA1998.jpg"" width="500">

### Data
The Sioux-Falls data used in the example is available from [repository](https://github.com/bstabler/TransportationNetworks).  You will need to unzip it and place in the "TrafficData" folder in order to run the code.  

## Technical Details
All relevand code is written in Julia using [JMP](https://github.com/JuliaOpt/JuMP.jl).  
* "fitTraffic.jl" and "trafficCVal.jl" perform most of work. The function "train" in "trafficCVal.jl" illustrates the key elements of fitting.
* The jupyter notebook "Traffic Experiments.ipynb" performs several of the validation experiments from the paper, but is sparsely documented.
* The file "genGraphs.r" creates all plots in the paper using R and ggplot2.  


