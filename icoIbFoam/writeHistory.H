if (runTime.outputTime())
{
    OFstream file
    (
        runTime.timePath()/"ib.dat"
    );

    const vectorField& ibC = ib.points();

    forAll(cylinderPoints, pointI)
    {
        label curPoint = cylinderPoints[pointI];

        scalar theta = ::atan2(ibC[curPoint].y(), ibC[curPoint].x());
        theta *= 180/M_PI;
        theta = 180 - theta;

        file << theta << tab << ibPressure[curPoint] << tab 
            << mag(ibSnGradU[curPoint]) << endl;
    }

    // Axis
    OFstream fileAxis
    (
        runTime.timePath()/"axis.dat"
    );

    const vectorField& C = mesh.Cf().internalField();
    vectorField Uf = fvc::interpolate(U)().internalField();

    forAll(axisLine, faceI)
    {
        label curFace = axisLine[faceI];

        fileAxis << C[curFace].x() << tab << Uf[curFace].x() << endl;
    }
}
