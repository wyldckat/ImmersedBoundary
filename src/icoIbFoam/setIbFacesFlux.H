{
    forAll (mesh.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            const immersedBoundaryFvPatch& ibPatch =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh.boundary()[patchI]
                );

            const labelList& ibFaces = ibPatch.ibFaces();

            forAll(ibFaces, faceI)
            {
                label curFace = ibFaces[faceI];

                if (curFace < mesh.nInternalFaces())
                {
                    phi.internalField()[curFace] = 
                        phiU.internalField()[curFace];
                }
                else
                {
                    label patchID = mesh.boundaryMesh().whichPatch(curFace);
                    label start = mesh.boundaryMesh()[patchID].start();
                    label localIndex = curFace - start;

                    phi.boundaryField()[patchID][localIndex] = 
                        phi.boundaryField()[patchID][localIndex];
                }
            }
        }
    }
}
