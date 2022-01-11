# lichen: Inverse modeling of lichen growth curves
Lichenometric dating represents a quick and affordable surface exposure dat-
ing method that has been widely used to provide a minimum age constraint
on tectonic and geomorphic landscape changes as well as buildings and an-
thropogenic landscape changes in various settings during the late Holocene.
Despite its widespread usage, this method has several limitations. Major
problems relate to the sampling of lichen population on any given rock sur-
face and the modeling of growth curves. In order to overcome these issues, it
has been suggested to subdivide the rock surface into some areas and mea-
sure the largest lichen thallus on each one. However, how to express the data
in terms of a probability distribution function and link it to an age of last
exposure of the rock surface are still a matter of debate. 
The code posted here is for the implementation of a novel approach to the 
modeling of lichen growth curves, which treats lichen growth as a continuous-time 
Markov process with a time-varying rate and additive Brownian noise. Given 
the growth rates, the probability distribution of the lichen population at 
any time can then be obtained by solving the Fokker-Planck equation. 
This method is illustrated using a dataset from the Huashan area of eastern 
China, which consists of measurements of the largest thalli on 12 rock surfaces 
of known age. We first build up the probability distribution of the lichen 
population for each rock surface based on extreme value theory and then use 
these to optimize the growth curve by minimizing the Jensen-Shannon divergence. 
A new method is also proposed to use the growth curve to map a sample of size 
data from an undated rock surface to the calendar age domain so as to yield 
a fully probabilistic estimate of the exposure age of the undated rock 
surface rather than a point estimate.
