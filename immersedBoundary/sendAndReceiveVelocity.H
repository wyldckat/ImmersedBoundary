    FieldField<Field, vector> procVelocities(Pstream::nProcs());
    forAll(procVelocities, procI)
    {
        procVelocities.set
        (
            procI,
            new vectorField(procCentres()[procI].size(), vector::zero)
        );
    }

    if (Pstream::parRun())
    {
//         const labelListList& procCells = *procCellsPtr_;

        // Send velocities
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField velocities(UI, procCells()[procI]);

                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI, 
                        velocities.size()*sizeof(vector)
                    );

                    toProc << velocities;
                }
            }
        }

        // Receive velocities
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI, 
                        procVelocities[procI].size()*sizeof(vector)
                    );

                    fromProc >> procVelocities[procI];
                }
            }
        }
    }
