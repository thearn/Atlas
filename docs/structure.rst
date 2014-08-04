============================================================
Structure of the Atlas plugin
============================================================

The Atlas plugin is organized in a hierarchical manner, separating
structural from aeronautical analysis through the use of OpenMDAO assemblies.

The `Aero()` assembly contains all components related to thrust, lift, drag, and
induced velocity calculations, which total 3 components.

The 'Structures()' assembly contains all components related to FEA, structural deformation,
and materials properties, which total 8 components.

The `AeroStructural()` assembly contains both a low-fidelity and high fidelity
`Aero()` instance, along with a `Structures()` instance. It also contains
a `DiscretizeProperties()` component, which takes discretization parameters
as input, and outputs variables discretized in space (typically along the
rotor blades).