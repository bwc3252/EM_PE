# Kilonova models

## Three-component analytic model

Currently the only kilonova model available is an analytic three-component model.
The model assumes each component is a homologously-expanding photosphere specified by an opacity, ejecta mass, ejecta velocity, and temperature floor.
The masses, velocities, and temperatures can all be fit, while the opacities are fixed.
See `Notes/notes.tex` for a more thorough description of the model.

This model is very fast due to being almost entirely analytic.
A typical PE run using this model takes about an hour if the likelihood is computed in 8 parallel processes on the machine used for testing.
