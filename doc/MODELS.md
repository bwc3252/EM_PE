# Kilonova models

## One-component analytic model

The simplest (and fastest) kilonova model.
This is an implementation of the model described in [Villar et al](https://arxiv.org/abs/1710.11576).
The model assumes a homologously-expanding photosphere specified by an opacity, ejecta mass, ejecta velocity, and temperature floor.
For simplicity the temperature floor is fixed to 4000 K, while the other parameters can vary.

This model is very fast due to being almost entirely analytic.

## Three-component analytic model

Based on the one-component analytic model, uses the sum of three independent components.

## Interpolated model
