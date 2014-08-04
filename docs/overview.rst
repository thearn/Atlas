============================================================
Overview of the Atlas model
============================================================

The Atlas plugin is an OpenMDAO implementation of the design code used to model
the AeroVelo Atlas quad-rotor human-powered helicopter. In 2012, the Atlas won the
Igor I. Sikorsky Human Powered Helicopter Competition by successfully hovering for
over 60 seconds at an altitude over 3 meters.

As a design problem, the objective is to minimize the total power required
to hover, with respect to the angular velocity of the rotors. Constraints relating
to structural dynamics (strain, buckling failure, material failure, etc) are
also included. The objective is optimized over two design points, one representing
conditions at low altitude (0.5 m) and high altitude (3.5 m).

Many additional design variables may also be activated and added to the problem,
such as rotor size, thickness of structural members, etc.