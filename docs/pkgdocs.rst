
================
Package Metadata
================

- **author:** Tristan A. Hearn

- **author-email:** tristan.a.hearn@nasa.gov

- **classifier**:: 

    Intended Audience :: Science/Research
    Topic :: Scientific/Engineering

- **description-file:** README.txt

- **entry_points**:: 

    [openmdao.component]
    Atlas.properties.DiscretizeProperties=Atlas.properties:DiscretizeProperties
    Atlas.properties.SparProperties=Atlas.properties:SparProperties
    Atlas.configuration.AtlasConfiguration=Atlas.configuration:AtlasConfiguration
    Atlas.thrust.ActuatorDiskInducedVelocity=Atlas.thrust:ActuatorDiskInducedVelocity
    Atlas.properties.QuadSparProperties=Atlas.properties:QuadSparProperties
    Atlas.aerostructural.AeroStructural=Atlas.aerostructural:AeroStructural
    Atlas.structures.Strains=Atlas.structures:Strains
    Atlas.lift_drag.LiftDrag=Atlas.lift_drag:LiftDrag
    Atlas.structures.FEM=Atlas.structures:FEM
    Atlas.aero.Aero2=Atlas.aero:Aero2
    Atlas.aero.Aero=Atlas.aero:Aero
    Atlas.properties.JointSparProperties=Atlas.properties:JointSparProperties
    Atlas.structures.Failures=Atlas.structures:Failures
    Atlas.properties.ChordProperties=Atlas.properties:ChordProperties
    Atlas.coefficients.DragCoefficient=Atlas.coefficients:DragCoefficient
    Atlas.structures.MassProperties=Atlas.structures:MassProperties
    Atlas.vortex.VortexRing=Atlas.vortex:VortexRing
    Atlas.thrust.Thrust=Atlas.thrust:Thrust
    Atlas.structures.Structures=Atlas.structures:Structures
    [openmdao.container]
    Atlas.properties.DiscretizeProperties=Atlas.properties:DiscretizeProperties
    Atlas.structures.Failure=Atlas.structures:Failure
    Atlas.properties.QuadSparProperties=Atlas.properties:QuadSparProperties
    Atlas.properties.JointSparProperties=Atlas.properties:JointSparProperties
    Atlas.lift_drag.Fblade=Atlas.lift_drag:Fblade
    Atlas.configuration.Flags=Atlas.configuration:Flags
    Atlas.structures.Strain=Atlas.structures:Strain
    Atlas.coefficients.DragCoefficient=Atlas.coefficients:DragCoefficient
    Atlas.vortex.VortexRing=Atlas.vortex:VortexRing
    Atlas.properties.SparProperties=Atlas.properties:SparProperties
    Atlas.configuration.AtlasConfiguration=Atlas.configuration:AtlasConfiguration
    Atlas.aerostructural.AeroStructural=Atlas.aerostructural:AeroStructural
    Atlas.structures.FEM=Atlas.structures:FEM
    Atlas.configuration.PrescribedLoad=Atlas.configuration:PrescribedLoad
    Atlas.thrust.Thrust=Atlas.thrust:Thrust
    Atlas.aero.Aero=Atlas.aero:Aero
    Atlas.properties.JointProperties=Atlas.properties:JointProperties
    Atlas.thrust.ActuatorDiskInducedVelocity=Atlas.thrust:ActuatorDiskInducedVelocity
    Atlas.structures.Structures=Atlas.structures:Structures
    Atlas.structures.BucklingFailure=Atlas.structures:BucklingFailure
    Atlas.lift_drag.LiftDrag=Atlas.lift_drag:LiftDrag
    Atlas.structures.MaterialFailure=Atlas.structures:MaterialFailure
    Atlas.structures.Strains=Atlas.structures:Strains
    Atlas.aero.Aero2=Atlas.aero:Aero2
    Atlas.structures.MassProperties=Atlas.structures:MassProperties
    Atlas.properties.ChordProperties=Atlas.properties:ChordProperties
    Atlas.structures.Failures=Atlas.structures:Failures

- **keywords:** openmdao

- **license:** Apache 2.0

- **maintainer:** Tristan A. Hearn

- **maintainer-email:** tristan.a.hearn@nasa.gov

- **name:** Atlas

- **requires-dist:** openmdao.main

- **requires-python**:: 

    >=2.6
    <3.0

- **static_path:** [ '_static' ]

- **summary:** OpenMDAO implementation of the Aerovelo Atlas human-powered helicopter design problem

- **version:** 0.1

