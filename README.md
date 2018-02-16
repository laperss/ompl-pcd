## PCD for OMPL
Probabilistic Cell Decomposition method implemented for the Open Motion-Planning Library (OMPL) [1]. 

PCD is a path planner made for high-dimensional static configuration spaces. [2]
Some code has been adapted from a previous implementation found [here](http://copp.cvs.sourceforge.net/). 

### Example
An example of a planned path: 

![Example path](https://user-images.githubusercontent.com/4593893/32837848-5a766340-ca0f-11e7-8747-9cda881bbd9e.png)

The computed cells are shown in red and green. The real obstacle is shown striped (black). A path found from this decomposition is shown in blue. 


![Example planning](https://people.kth.se/~laperss/assets/images/animation_blocks_example.gif)



## References
[1] http://ompl.kavrakilab.org/

[2] http://ieeexplore.ieee.org/document/1307193/
