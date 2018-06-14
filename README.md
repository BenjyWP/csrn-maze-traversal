[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)

# Maze Traversal with Extended Kalman Filter

Cellular simultaneous recurrent neural network (SRN) has been shown to be a function approximator more powerful than the multilayer perceptron (MLP). This means that the complexity of MLP would be prohibitively large for some problems while SRN could realize the desired mapping with acceptable computational constraints. The speed of training of complex recurrent networks is crucial to their successful application. This work improves the previous results by training the network with extended Kalman filter (EKF). We implemented a generic cellular SRN (CSRN) and applied it for solving two challenging problems: 2-D maze navigation and a subset of the connectedness problem. The speed of convergence has been improved by several orders of magnitude in comparison with the earlier results in the case of maze navigation, and superior generalization has been demonstrated in the case of connectedness.

The implications of this improvement are discussed in [1].

The above work has been extended to perform image affine transformation using CSRN architecture. The details of this work can be found in [2]. Moreover, an improved learning algorithm known as unscented Kalman filter (UKF) has been proposed in [3].

## Publications

The following publications cover different stages of the development of this project.

1. R. Ilin, R. Kozma, and P. Werbos, "Beyond Feedforward Models Trained by Backpropogation: a Practical Tool for a More Efficient Universal Approximator," in Proc. IEEE Neural Networks Transactions on, May 2008

2. R. Ilin, R. Kozma, and P. Werbos, “Cellular SRN trained by extended Kalman filter shows promise for ADP,” in Proc. IEEE IJCNN, Jul. 2006, pp. 506–510.

3. J. K. Anderson and K. M. Iftekharuddin, “Learning topological image transforms using cellular simultaneous recurrent networks,” Neural Networks (IJCNN), The 2013 International Joint Conference on, 2013, pp. 1-9.

4. L. Vidyaratne, M. Alam, J. K. Anderson, and K. M. Iftekharuddin, “Improved training of cellular SRN using Unscented Kaiman Filtering for ADP,” in Neural Networks (IJCNN), 2014 International Joint Conference on, 2014, pp. 993-1000.

## License

This project is licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License](https://creativecommons.org/licenses/by-nc-sa/4.0/) - see the [LICENSE.md](LICENSE.md) file for details.
