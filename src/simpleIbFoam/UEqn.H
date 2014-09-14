    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    UEqn().relax();

    eqnResidual = solve
    (
        UEqn() == -cellIbMask*fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
