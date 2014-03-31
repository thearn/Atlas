
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
    Atlas.properties.SparProperties=Atlas.properties:SparProperties
    Atlas.properties.DiscretizeProperties=Atlas.properties:DiscretizeProperties
    Atlas.lift_drag.LiftDrag=Atlas.lift_drag:LiftDrag
    Atlas.thrust.ActuatorDiskInducedVelocity=Atlas.thrust:ActuatorDiskInducedVelocity
    Atlas.thrust.Thrust=Atlas.thrust:Thrust
    Atlas.vortex.InducedVelocity=Atlas.vortex:InducedVelocity
    Atlas.structures.Strains=Atlas.structures:Strains
    Atlas.aero.Configuration=Atlas.aero:Configuration
    Atlas.aero.AeroCalc=Atlas.aero:AeroCalc
    Atlas.structures.MassProperties=Atlas.structures:MassProperties
    Atlas.properties.ChordProperties=Atlas.properties:ChordProperties
    Atlas.coefficients.DragCoefficient=Atlas.coefficients:DragCoefficient
    Atlas.structures.Failures=Atlas.structures:Failures
    Atlas.vortex.VortexRing=Atlas.vortex:VortexRing
    Atlas.structures.FEM=Atlas.structures:FEM
    Atlas.structures.Structures=Atlas.structures:Structures
    [openmdao.container]
    Atlas.properties.DiscretizeProperties=Atlas.properties:DiscretizeProperties
    Atlas.aero.Configuration=Atlas.aero:Configuration
    Atlas.structures.Failure=Atlas.structures:Failure
    Atlas.structures.PrescribedLoad=Atlas.structures:PrescribedLoad
    Atlas.lift_drag.Fblade=Atlas.lift_drag:Fblade
    Atlas.structures.Strain=Atlas.structures:Strain
    Atlas.structures.JointProperties=Atlas.structures:JointProperties
    Atlas.vortex.VortexRing=Atlas.vortex:VortexRing
    Atlas.coefficients.DragCoefficient=Atlas.coefficients:DragCoefficient
    Atlas.properties.SparProperties=Atlas.properties:SparProperties
    Atlas.vortex.InducedVelocity=Atlas.vortex:InducedVelocity
    Atlas.structures.FEM=Atlas.structures:FEM
    Atlas.thrust.Thrust=Atlas.thrust:Thrust
    Atlas.aero.AeroCalc=Atlas.aero:AeroCalc
    Atlas.structures.Flags=Atlas.structures:Flags
    Atlas.thrust.ActuatorDiskInducedVelocity=Atlas.thrust:ActuatorDiskInducedVelocity
    Atlas.structures.Structures=Atlas.structures:Structures
    Atlas.structures.BucklingFailure=Atlas.structures:BucklingFailure
    Atlas.lift_drag.LiftDrag=Atlas.lift_drag:LiftDrag
    Atlas.structures.MaterialFailure=Atlas.structures:MaterialFailure
    Atlas.structures.Strains=Atlas.structures:Strains
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

