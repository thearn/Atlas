
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
    Atlas.vortex.inducedVelocity=Atlas.vortex:inducedVelocity
    Atlas.failure.failureCalc=Atlas.failure:failureCalc
    Atlas.properties.WireProperties=Atlas.properties:WireProperties
    Atlas.structures.Structural=Atlas.structures:Structural
    Atlas.structures.MassProperties=Atlas.structures:MassProperties
    Atlas.coefficients.dragCoefficientFit=Atlas.coefficients:dragCoefficientFit
    Atlas.properties.ChordProperties=Atlas.properties:ChordProperties
    Atlas.coefficients.dragCoefficient=Atlas.coefficients:dragCoefficient
    Atlas.failure.torsionalBucklingFailure=Atlas.failure:torsionalBucklingFailure
    Atlas.vortex.vortexRing=Atlas.vortex:vortexRing
    Atlas.structures.FEM=Atlas.structures:FEM
    [openmdao.container]
    Atlas.properties.DiscretizeProperties=Atlas.properties:DiscretizeProperties
    Atlas.properties.SparProperties=Atlas.properties:SparProperties
    Atlas.structures.Flags=Atlas.structures:Flags
    Atlas.structures.PrescribedLoad=Atlas.structures:PrescribedLoad
    Atlas.vortex.inducedVelocity=Atlas.vortex:inducedVelocity
    Atlas.failure.failureCalc=Atlas.failure:failureCalc
    Atlas.properties.WireProperties=Atlas.properties:WireProperties
    Atlas.structures.FBlade=Atlas.structures:FBlade
    Atlas.coefficients.dragCoefficient=Atlas.coefficients:dragCoefficient
    Atlas.structures.Structural=Atlas.structures:Structural
    Atlas.structures.FEM=Atlas.structures:FEM
    Atlas.structures.MassProperties=Atlas.structures:MassProperties
    Atlas.coefficients.dragCoefficientFit=Atlas.coefficients:dragCoefficientFit
    Atlas.properties.ChordProperties=Atlas.properties:ChordProperties
    Atlas.structures.JointProperties=Atlas.structures:JointProperties
    Atlas.vortex.vortexRing=Atlas.vortex:vortexRing
    Atlas.failure.torsionalBucklingFailure=Atlas.failure:torsionalBucklingFailure

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

