    FieldField<Field, Type> procPsi(Pstream::nProcs());
    forAll(procPsi, procI)
    {
        procPsi.set
        (
            procI,
            new Field<Type>(procCentres()[procI].size(), pTraits<Type>::zero)
        );
    }

    if (Pstream::parRun())
    {
        // Send
        for (label procI=0; procI<Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                Field<Type> psi_(psiI, procCells()[procI]);
                
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI, 
                        psi_.size()*sizeof(Type)
                    );

                    toProc << psi_;
                }
            }
        }

        // Receive
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
                        procPsi[procI].size()*sizeof(Type)
                    );

                    fromProc >> procPsi[procI];
                }
            }
        }
    }
