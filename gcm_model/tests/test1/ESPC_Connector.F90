!
!
#define FILENAME "ESPC_Connector.F90"
#define CONTEXT  line=__LINE__,file=__FILE__


module  ESPC_Connector

  use ESMF
  use NUOPC

  use NUOPC_Connector, only: &
    con_routine_SS      => SetServices, &
    con_label_ExecuteRH => label_ExecuteRouteHandle, &
    con_label_ReleaseRH => label_ReleaseRouteHandle, &
    con_label_ComputeRH => label_ComputeRouteHandle, &
    NUOPC_ConnectorGet, NUOPC_ConnectorSet

  implicit none

  private

  public SetServices

  integer nPets,lPet

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(connector, rc)
    type(ESMF_CplComp)  :: connector
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! the NUOPC connector component will register the generic methods
    call NUOPC_CompDerive(connector, con_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, rcToReturn=rc)) RETURN

    ! attach specializing method to compute the connection RouteHandle
    call NUOPC_CompSpecialize(connector, specLabel=con_label_ComputeRH, &
      specRoutine=ComputeRH, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, rcToReturn=rc)) RETURN

    call NUOPC_CompSpecialize(connector, specLabel=con_label_ExecuteRH, &
      specRoutine=ExecuteRH, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, rcToReturn=rc)) RETURN

    call NUOPC_CompSpecialize(connector, specLabel=con_label_ReleaseRH, &
      specRoutine=ReleaseRH, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, rcToReturn=rc)) RETURN

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine ComputeRH(connector, rc)
    type(ESMF_CplComp)  :: connector
    CHARACTER(len=10) :: cplName
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call ESMF_CplCompGet(connector,petCount=nPets,localPet=lPet,name=cplName,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, rcToReturn=rc)) RETURN


    if(lPet.eq.0) print *,"ComputeRH called... ",cplName

  end subroutine ComputeRH

  !-----------------------------------------------------------------------------

  subroutine ExecuteRH(connector, rc)
    type(ESMF_CplComp)  :: connector
    CHARACTER(len=10) :: cplName
    integer, intent(out) :: rc

    call ESMF_CplCompGet(connector,name=cplName,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, rcToReturn=rc)) RETURN

    if(lPet.eq.0) print *,"ExecuteRH called... ", cplName

  end subroutine ExecuteRH

  subroutine ReleaseRH(connector, rc)
    type(ESMF_CplComp)  :: connector
    integer, intent(out) :: rc
    CHARACTER(len=10) :: cplName

    call ESMF_CplCompGet(connector,name=cplName,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, CONTEXT, rcToReturn=rc)) RETURN

    if(lPet.eq.0) print *,"ReleaseRH called ",cplName

  end subroutine ReleaseRH

!===============================================
!========================================================
!================================================================

end module

