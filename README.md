# IBIOColorDetect

An isetbio computational observer for threshold detection thresholds of color Gabor stimuli.  The goal of these calculations is to understand how various factors early in the visual pathways limit performance on this simple and rigorously characterized visual task.

The place to get started is with the tutorials.  These are designed to build up to the full caclulation in stages.  Code that is illustrated in the initial tutorials is encapsulated in routines in the toolbox, with those routines called in more advanced tutorials.

## Tutorial List

The tutorials live in the tutorials folder of the respository.  They are designed to provide an introduction to the ideas and code.

t_colorGaborScene - Shows how to make an isetbio scene representing a colored Gabor pattern presented on a calibrated CRT monitor.  This is the basic stimulus whose detection threshold we are modeling in this project.  The code illustrated in this tutorial is encapsulated in function colorGaborSceneCreate.

t_colorGaborConeAbsorptionMovie - Shows how to take a temporally windowed color Gabor stimulus (Gaussian window) and compute a movie of the cone mosaic isomerizations at each time sampling point.  This relies on function colorGaborSceneCreate.

t_colorGaborConeCurrentEyeMovementsMovie - This goes further and adds modeling of eye movements as well as the tranformation from absorption to the outer segment cone current.  The functionality showin in this tutorial is encapsulated in the function colorDetectResponseInstanceConstruct.

t_colorGaborConeCurrentEyeMovementsResponseInstances - For building classifiers, we need to get multiple noisy instances of the responses to stimuli.  This shows how to get such instances, by calling function colorDetectResponseInstanceConstruct.  It does so for multiple contrasts and color directions, and saves the instances for use by the colorGaborDetectFindPerformance tutorial immediately below.  The output can be rather large.  See comment at top of file for a brief description of where the output goes and how to control it.

t_colorGaborDetectFindPerformance - Use an SVM as a computational observer to find thresholds.  This reads the output of t_colorGaborConeCurrentEyeMovementsResponseInstances, and saves its output for t_plotColorGaborDetectThresholdsOnLMPlane.

t_colorGaborDetectIllustrateClassificationBoundary - Illustrate the classification boundary for a specified case of the output of t_colorDetectFindPerformance.ß

t_plotColorGaborDetectThresholdsOnLMPlane - Fit psychometric function to the output of t_colorGaborDetectFindPerformance and use this to plot an isodetection threshold contour in the LM plane, with an ellipse fit to it.  This could be pretty easily generalized to other planes or full ellipsoids.

t_colorThresholdEllipsoids - This tutorial implements the model developed in Poirson & Wandell (1996) and allows one to visual color threshold ellipsoids for various choices of spatial frequency, according to that model.  These threshold ellipsoids summarize one of the data sets we would like to model with the code in this repository.

t_colorThresholdEllipsoidFit - Illustrates how to fit an ellipsoid to color threshold data for a variety of directions in color space.
