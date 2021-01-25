## Modules

The 3 main modules, **Score Peaks**, **Find and Score Peaks**, and **Find and Score Peaks with Input**, are found here, along with the LanceOtron.py script of classes and functions used, model weights for the neural network, and coefficients for the standard scaler normilisation.

### TensorFlow versions

The scripts here use TensorFlow 2 using Keras as a high-level API, implemented as *module of TensorFlow*. While it is possible to save the model structure and weights as a single unit, backwards compatibility becomes and issue. For this reason we chose to include a 'build_model' function in the scripts, which also explicity details the model structure, and save only the model weights. To make this compatible with TensorFlow 1, the Keras python package would need to be loaded and aliased appropriately, but could be accomplished relatively easily.
