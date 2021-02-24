MODULE nrtype
    !====================================================================
    !Define the types used for the program
    !====================================================================
    INTEGER, PARAMETER :: SP = KIND(1.0)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: SPC = KIND((1.0_sp,1.0_sp))
    INTEGER, PARAMETER :: DPC = KIND((1.0_dp,1.0_dp))
    INTEGER, PARAMETER :: LGT = KIND(.true.)

END MODULE nrtype

