    FieldField<Field, scalar> procPressures(Pstream::nProcs());
    forAll(procPressures, procI)
    {
        procPressures.set
        (
            procI,
            new scalarField(procCentres()[procI].size(), 0)
        );
    }

    if (Pstream::parRun())
    {
//         const labelListList& procCells = *procCellsPtr_;

        // Send pressure
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                scalarField pressures(pI, procCells()[procI]);
                
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI, 
                        pressures.size()*sizeof(scalar)
                    );

                    toProc << pressures;
                }
            }
        }

        // Receive pressures
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
                        procPressures[procI].size()*sizeof(vector)
                    );

                    fromProc >> procPressures[procI];
                }
            }
        }
    }
