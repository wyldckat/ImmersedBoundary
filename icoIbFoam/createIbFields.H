vectorField ibVelocity(ib.points().size(), vector::zero);
scalarField ibSnGradP(ib.points().size(), 0.0);

scalarField ibPressure(ib.points().size(), 0);
vectorField ibSnGradU(ib.points().size(), vector::zero);
