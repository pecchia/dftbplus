!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

module phonons_initphonons
  use dftbp_common_accuracy
  use dftbp_common_atomicmass
  use dftbp_common_constants
  use dftbp_common_environment
  use dftbp_common_globalenv
  use dftbp_common_unitconversion
  use dftbp_dftb_periodic
  use dftbp_extlibs_negf, only : z_CSR, create, destroy
  use dftbp_io_charmanip
  use dftbp_io_fileid
  use dftbp_io_hsdparser, only : parseHSD, dumpHSD
  use dftbp_io_hsdutils
  use dftbp_io_hsdutils2
  use dftbp_io_message
  use dftbp_io_tokenreader
  use dftbp_io_xmlutils
  use dftbp_math_simplealgebra
  use dftbp_transport_negfvars
  use dftbp_type_linkedlist
  use dftbp_type_oldskdata
  use dftbp_type_typegeometryhsd
  use dftbp_type_wrappedintr
  use xmlf90_flib_dom
  use phonons_mesh
  implicit none 
  private

  character(len=*), parameter :: rootTag = "phonons"
  character(len=*), parameter :: autotestTag = "autotest.tag"
  character(len=*), parameter :: hsdInput = "phonons_in.hsd"
  character(len=*), parameter :: hsdParsedInput = "phonons_pin.hsd"
  character(len=*), parameter :: xmlInput = "phonons_in.xml"
  character(len=*), parameter :: xmlParsedInput = "phonons_pin.xml"

  public :: initProgramVariables, destructProgramVariables
  public :: TPdos, autotestTag

  type TPdos
    type(TWrappedInt1), allocatable :: iAtInRegion(:)
    character(lc), allocatable :: regionLabels(:)
  end type TPdos
  
  !> Identity of the run
  integer, public :: identity

  !> Geometry container
  type(TGeometry), public :: geo

  !> Projected dos infos
  type(TPdos), public :: pdos

  !> Container of transport parameters
  type(TTransPar), public :: transpar

  !> Container of tunneling parameters
  type(TNEGFtundos), public :: tundos

  !> verbose flag
  logical, public :: tVerbose

  !> expose z_CSR type used here
  public :: z_CSR
  !> Core Variables 
  real(dp), allocatable, public :: atomicMasses(:)
  real(dp), allocatable, public :: atomicMassesPert(:)
  real(dp), allocatable, public :: dynMatrix(:,:)
  type(z_CSR), public, target :: dynMatCsr
  real(dp), allocatable, public :: dynMatrixPert(:,:)
  type(z_CSR), public, target :: dynMatPertCsr
  integer, allocatable, public :: iMovedAtoms(:)
  integer, allocatable, public :: centralCellAtom(:)
  integer, public :: nMovedAtom, nAtomUnitCell
  !> Maps of primitive <-> supercell atom indices
  type(TWrappedInt1), allocatable, public :: p2s_map(:)    
  integer, allocatable, public :: s2p_map(:)

  !> Kpoints information
  real(dp), allocatable, public :: kPoint(:,:), kWeight(:)
  integer, public :: nKPoints
  type(TMesh), public :: kmesh


  !> maps atom index in central cell
  integer, allocatable, public  :: Img2CentCell(:)

  !> neighbor list
  type(TNeighbourList), public :: neighbourList

  !> number of neighbors
  integer, allocatable, target, public :: nNeighbour(:)

  !> cutoff for Hessian
  real(dp), public :: cutoff

  !> Temperature range
  real(dp), public :: TempMin, TempMax, TempStep

  !> Modes to analyze (e.g., longitudinal, transverse, in-plane, etc)
  integer, public :: selTypeModes

  !> order=2 means harmonic, 3 is anharmonic 3rd order, etc.
  integer, public :: order

  !> Atomic temperature
  real(dp), public :: atTemperature

  !> Whether modes should be computed
  logical, public :: tCompModes

  !> Whether modes should be plotted
  logical, public :: tPlotModes

  !> Whether modes should be animated
  logical, public :: tAnimateModes

  !> creates output for xmakemol
  logical, public :: tXmakeMol

  !> whether transport calculations should be done 
  logical, public :: tTransport

  !> whether phonon dispersions should be computed
  logical, public :: tPhonDispersion

  !> whether phonon density of states should be computed
  logical, public :: tDensityOfStates

  !> whether a kmesh has been defined 
  logical, public :: tKMesh

  !> compute phonon relaxation times with T-matrix
  logical, public :: tRelaxationTimes

  !> number of repeated cells along lattice vectors
  integer, public :: nCells(3)

  !> volume or Area (2D) of the supercell geometry
  real(dp), public :: cellVolume

  !> whether phonon dispersions should be computed
  character(4), public :: outputUnits

  !> whether taggedoutput should be written
  logical, public :: tWriteTagged

  !> Which phonon modes to animate
  integer, allocatable, public :: modesToPlot(:)

  !> Number of phonon modes to animate
  integer, public :: nModesToPlot

  !> number of cycles in mode animation
  integer, public :: nCycles

  !> number of steps in mode animation
  integer, public, parameter :: nSteps = 10

  !> Version of the current parser
  integer, parameter :: parserVersion = 4

  !> Version of the oldest parser, for which compatibility is maintained
  integer, parameter :: minVersion = 4

  !> Delta broadening for Green's functions
  real(dp), public :: delta

  !> Container type for parser related flags.
  type TParserFlags

    !> stop after parsing?
    logical :: tStop

    !> Continue despite unprocessed nodes
    logical :: tIgnoreUnprocessed

    !> write parsed HSD input
    logical :: tWriteHSD

    !> write TaggedOutput

    logical :: tWriteTagged
  end type TParserFlags


  !> constants parameters
  type TModeEnum
    !> along x,y,z
    integer :: ALLMODES = 0
    !> along x
    integer :: XX = 1
    !> along y
    integer :: YY = 2
    !> along z
    integer :: ZZ = 3
    !> along z (assume transport along z)
    integer :: LONGITUDINAL = 4
    !> along x,y
    integer :: TRANSVERSE = 5
    !> along x,z (assume 2D structure in x-z)
    integer :: INPLANE = 6
    !> along y
    integer :: OUTOFPLANE = 7
  end type TModeEnum

  type(TModeEnum), public, parameter :: modeEnum = TModeEnum()

contains

  !> Initialise program variables
  subroutine initProgramVariables(env)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    ! locals
    type(fnode), pointer :: hsdTree, root, node, tmp
    type(fnode), pointer :: child, child2, value
    type(fnodeList), pointer :: children
    type(string) :: buffer, buffer2, modif
    integer :: inputVersion
    integer :: ii, jj, iSp1, iAt
    logical :: tHSD, reqMass, tBadKPoints
    real(dp), allocatable :: speciesMass(:)
    integer, allocatable :: tmpI1(:)
    real(dp) :: drop, rTmp
    integer :: nDerivs, nGroups
    type(TParserflags) :: parserFlags
    type(TListIntR1) :: li1

    integer :: cubicType, quarticType
    integer :: idxCell(2), nn
    integer :: ip, is

    write(stdOut, "(/, A)") "Starting initialization..."
    write(stdOut, "(A80)") repeat("-", 80)

    nGroups = 1
#:if WITH_MPI
    call env%initMpi(nGroups)
#:endif

    !! Read in input file as HSD or XML.
    write(stdOut, "(A)") "Interpreting input file '" // hsdInput // "'"
    write(stdOut, "(A)") repeat("-", 80)
    call parseHSD(rootTag, hsdInput, hsdTree)
    call getChild(hsdTree, rootTag, root)

    !! Check if input version is the one, which we can handle
    !! Handle parser options
    call getChildValue(root, "Options", tmp, "", child=child, &
        &list=.true., allowEmptyValue=.true.)
    call readOptions(child, root, parserFlags)

    call getChild(root, "Geometry", tmp)
    call readGeometry(tmp, geo)

    ! Read Transport block
    ! This defines system partitioning
    call getChild(root, "Transport", child, requested=.false.)
    if (associated(child)) then
      tTransport = .true.
      call readTransportGeometry(child, geo, transpar)
    else
      tTransport = .false.
    end if

    call getChildValue(root, "Atoms", buffer2, "1:-1", child=child)
    call getSelectedAtomIndices(child, char(buffer2), geo%speciesNames, geo%species, iMovedAtoms)
    nMovedAtom = size(iMovedAtoms)

    call getChild(root, "ComputeModes",child=node,requested=.false.)
    if (associated(node)) then
      tCompModes = .true.
      call getChild(root, "DisplayModes",child=node,requested=.false.)
      if (associated(node)) then
        tPlotModes = .true.
        call getChildValue(node, "PlotModes", buffer2, "1:-1", child=child, &
            &multiple=.true.)
        call getSelectedIndices(child, char(buffer2), [1, 3 * nMovedAtom], modesToPlot)
        nModesToPlot = size(modesToPlot)
        call getChildValue(node, "Animate", tAnimateModes, .true.)
        call getChildValue(node, "XMakeMol", tXmakeMol, .true.)
      else
        nModesToPlot = 0
        tPlotModes = .false.
        tAnimateModes = .false.
        tXmakeMol = .false.
      end if
    else
      tCompModes = .false.
    end if

    if (tAnimateModes.and.tXmakeMol) then
      nCycles = 1
    else
      nCycles = 3
    end if

    ! Reading K-points for Phonon Dispersion calculation
    call getChild(root, "PhononDispersion", child=node, requested=.false.)
    if  (associated(node))  then
      tPhonDispersion = .true.
      print*,'Phonon bandstructure calculation'
      if (.not.geo%tPeriodic) then
         call detailedError(root,"Bandstructure calculation with non-periodic structure")
      end if
      call init(li1)
      call getChildValue(node, "supercell", 3, li1)
      call asVector(li1, nCells)
      call destruct(li1)
      nAtomUnitCell = geo%nAtom/(nCells(1)*nCells(2)*nCells(3))
      write(stdOut, *) 'Number of atoms in Unit Cell:', nAtomUnitCell
      call getChildValue(node, "centralCellAtoms", buffer, child=child2, multiple=.true.)

      print*, char(buffer)

      if (associated(child2)) then    
        call getSelectedAtomIndices(child2, char(buffer), geo%speciesNames, geo%species, &
              & centralCellAtom)
        if (size(centralCellAtom) /= nAtomUnitCell) then
           call detailedError(node,"centralCellAtoms inconsistent with supercell definition")
        end if
        if (any(centralCellAtom<1) .or. any(centralCellAtom>geo%nAtom)) then
           call detailedError(node,"centralCellAtoms outside the structure range")
        end if
        call getMaps(centralCellAtom, nCells, p2s_map, s2p_map)
        print*, 'primitive -> supercell'
        do ip = 1, size(centralCellAtom)
           print*,ip,':',p2s_map(ip)%data
        end do
        print*, 'supercell -> primitive'
        do is = 1, geo%nAtom
          print*, is, s2p_map(is)
        end do
      end if
      tAnimateModes = .false.
      tXmakeMol = .false.
      call getChild(node, "DisplayModes",child=child,requested=.false.)
      if (associated(child)) then
        tPlotModes = .true.
        call getChildValue(child, "PlotModes", buffer2, "1:-1", child=child, &
            &multiple=.true.)
        call getSelectedIndices(child, char(buffer2), [1, 3 * nMovedAtom], modesToPlot)
        nModesToPlot = size(modesToPlot)
        ! checks that modesToPlot is consitent with unit cell
        if (nModesToPlot > 3*nAtomUnitCell .or. any(modesToPlot .gt. 3*nAtomUnitCell)) then
          write(stdOut,*) 'Maximum expected number of modes:',3*nAtomUnitCell    
          call detailedError(node,"PlotModes inconsistent with atoms in unit cell")
        end if
      else
        nModesToPlot = 0
        tPlotModes = .false.
      end if
      print*,'Plot Modes:',nModesToPlot
      call getChildValue(node, "outputUnits", buffer, "H")
      select case(trim(char(buffer)))
      case("H", "eV" , "meV", "THz", "cm")
        outputUnits=trim(char(buffer))
      case default
        call detailedError(node,"Unknown outputUnits "//trim(char(buffer)))
      end select
      print*,'Read k-points:'
      call readKPoints(node, geo, tKMesh, tBadKpoints)
    else
      tPhonDispersion = .false.
    end if
    

    ! Read the atomic masses from SlaterKosterFiles or Masses
    allocate(speciesMass(geo%nSpecies))
    call getChild(root,"Masses",child=node, requested=.true.)
    if ( associated(node) ) then
      call getChild(node, "SlaterKosterFiles", child=value,requested=.false.)
      if ( associated(value) ) then
        call readSKfiles(value, geo, speciesMass)
      else
        call readMasses(node, geo, speciesMass)
      endif
    endif
    allocate(atomicMasses(nMovedAtom))
    do iAt = 1, nMovedAtom
      atomicMasses(iAt) = speciesMass(geo%species(iMovedAtoms(iAt)))
    end do
    deallocate(speciesMass)

    ! --------------------------------------------------------------------------------------
    ! Reading Hessian block parameters
    ! --------------------------------------------------------------------------------------
    call getChild(root, "Hessian", child=node, requested=.true.)
    ! cutoff used to cut out interactions
    call getChildValue(node, "Cutoff", cutoff, 9.45_dp, modifier=modif, child=value)
    call convertByMul(char(modif), lengthUnits, value, cutoff)
    ! Sets directly a drop off in the Hessian/DynMat matrix elements
    call getChildValue(node, "Drop", drop, 1.d-30, modifier=modif, child=value) 

    ! Reading the actual Hessian matrix
    call getChildValue(node, "Matrix", value, child=child)
    call getNodeName(value, buffer)
    select case(trim(char(buffer)))
    case ("dftb")
      call readDftbHessian(value, atomicMasses, DynMatrix)
    case ("cp2k")
      call readCp2kHessian(value, atomicMasses, DynMatrix)
    case ("vasp")
      call readVaspHessian(value, atomicMasses, DynMatrix)
    case ("dynmatrix")
      call readDynMatrix(value, DynMatrix)
    case ("lammps")
      call readLammpsMatrix(value, DynMatrix)
    case default
      call detailedError(node,"Unknown Hessian type "//char(buffer))
    end select

    ! CSR transformation for interfacing with libNEGF
    if (tTransport) then
      call dns2csr(dynMatrix, dynMatCsr, drop)
      deallocate(dynMatrix)
    end if

    call getChildValue(node, "DensityOfStates", tDensityOfStates, .false.) 
    if (tDensityOfStates .and. .not.tKMesh) then
      call detailedError(node,"DensityOfStates requires definition of a kmesh")     
    end if
    
    call getChildValue(node, "PhononLifetimes", tRelaxationTimes, .false.)
    if (tRelaxationTimes .and. .not.tKMesh) then
      call detailedError(node,"PhononLifetimes require definition of a kmesh")     
    end if
    if (tRelaxationTimes .and. .not.tPhonDispersion .and. .not.tPlotModes) then
      call detailedError(node,"Phonon lifetimes possible only with PhononDispersion{} &
            & and DisplayModes{}")  
    end if

    if (tRelaxationTimes) then
      atomicMassesPert = atomicMasses
      call getChildren(node, "PerturbedMasses", children, .false.)
      do ii = 1, getLength(children)
        call getItem1(children, ii, child)
        call getChildValue(child, "Atoms", buffer, child=child2, multiple=.true.)
        call getSelectedAtomIndices(child2, char(buffer), geo%speciesNames, geo%species, tmpI1)
        call getChildValue(child, "Mass", rTmp)
        do jj = 1, size(tmpI1)
          atomicMassesPert(tmpI1(jj)) = rTmp * amu__au
        end do
        deallocate(tmpI1) 
      end do

      call getChildValue(node, "MatrixPerturbed", value, child=child)
      call getNodeName(value, buffer)
      select case(trim(char(buffer)))
      case ("dftb")
        call readDftbHessian(value, atomicMassesPert, DynMatrixPert)
      case ("cp2k")
        call readCp2kHessian(value, atomicMassesPert, DynMatrixPert)
      case ("vasp")
        call readVaspHessian(value, atomicMassesPert, DynMatrixPert)
      case ("dynmatrix")
        call readDynMatrix(value, DynMatrixPert)
      case ("lammps")
        call readLammpsMatrix(value, DynMatrixPert)
      case default
        call detailedError(node,"Unkown Hessian type "//char(buffer))
      end select
      call getChildValue(node, "Delta", tundos%delta, &
          &0.0001_dp, modifier=modif, child=value)
      call convertByMul(char(modif), energyUnits, value, tundos%delta)
      delta = tundos%delta
    end if       

    ! --------------------------------------------------------------------------------------
    ! Reading cubic forces
    ! --------------------------------------------------------------------------------------
    order = 2
    call getChild(root, "Cubic", child=child, requested=.false.)
    if (associated(child)) then
      order = 3
      call getChildValue(child, "Matrix", buffer, child=child)
      select case(trim(char(buffer)))
      case ("gaussian")
        cubicType = 1
      case default
        call detailedError(root,"Unknown Cubic forces type "//char(buffer))
      end select
    end if

    call buildNeighbourList()

    call getChildValue(root, "Analysis", tmp, "", child=child, list=.true., &
        &allowEmptyValue=.true., dummyValue=.true.)

    if (associated(tmp)) then
      if (tPhonDispersion) then
         call detailedError(root, "Analysis and PhononDispersion cannot coexist")
      end if
      call readAnalysis(child, geo, pdos, tundos, transpar, atTemperature)
    endif

    !! Issue warning about unprocessed nodes
    write(stdOut, "(/, A)") "check unprocessed nodes..."
    call warnUnprocessedNodes(root, parserFlags%tIgnoreUnprocessed )

    !! Dump processed tree in HSD and XML format
    if (tIoProc .and. parserFlags%tWriteHSD) then
      call dumpHSD(hsdTree, hsdParsedInput)
      write(stdOut, '(/,/,A)') "Processed input in HSD format written to '" &
          &// hsdParsedInput // "'"
    end if

    !! Stop, if only parsing is required
    if (parserFlags%tStop) then
      call error("Keyword 'StopAfterParsing' is set to Yes. Stopping.")
    end if

    write(stdOut, "(/, A)") "Initialization done..."

  end subroutine initProgramVariables

  !!* destruct the program variables created in initProgramVariables
  subroutine destructProgramVariables()
    deallocate(atomicMasses)
    if (allocated(dynMatrix)) then
      deallocate(dynMatrix)
    end if  
    if (allocated(dynMatrixPert)) then
      deallocate(dynMatrixPert)
    end if 
    if (allocated(dynMatCsr%nzval)) then 
      call destroy(dynMatCsr)
    end if  
    if (allocated(iMovedAtoms)) then
      deallocate(iMovedAtoms)
    end if
    if (allocated(modesToPlot)) then
      deallocate(modesToPlot)
    end if
    write(stdOut, "(/,A)") repeat("=", 80)
  end subroutine destructProgramVariables


  !!* Read in parser options (options not passed to the main code)
  !!* @param node Node to get the information from
  !!* @param root Root of the entire tree (in the case it must be converted)
  !!* @param flags Contains parser flags on exit.
  subroutine readOptions(node, root, flags)
    type(fnode), pointer :: node
    type(fnode), pointer :: root
    type(TParserFlags), intent(out) :: flags

    integer :: inputVersion
    type(fnode), pointer :: child

    !! Check if input needs compatibility conversion.
    call getChildValue(node, "ParserVersion", inputVersion, parserVersion, &
        &child=child)
    if (inputVersion < 1 .or. inputVersion > parserVersion) then
      call detailedError(child, "Invalid parser version (" // i2c(inputVersion)&
          &// ")")
    elseif (inputVersion < minVersion) then
      call detailedError(child, &
          &"Sorry, no compatibility mode for parser version " &
          &// i2c(inputVersion) // " (too old)")
    end if

    call getChildValue(node, "WriteAutotestTag", tWriteTagged, .false.)
    call getChildValue(node, "WriteHSDInput", flags%tWriteHSD, .true.)
    call getChildValue(node, "StopAfterParsing", flags%tStop, .false.)

    call getChildValue(node, "IgnoreUnprocessedNodes", &
        &flags%tIgnoreUnprocessed, .false.)

  end subroutine readOptions

  !!* Read in the geometry stored as xml in internal or gen format.
  !!* @param geonode Node containing the geometry
  !!* @param geo     Contains the geometry information on exit
  subroutine readGeometry(geonode, geo)
    type(fnode), pointer :: geonode
    type(TGeometry), intent(out) :: geo

    type(fnode), pointer :: child, value
    type(string) :: buffer

    call getChildValue(geonode, "", value, child=child)
    call getNodeName(value, buffer)
    select case (char(buffer))
    case ("genformat")
      call readTGeometryGen(value, geo)
    case default
      call setUnprocessed(value)
      call readTGeometryHSD(child, geo)
    end select

    if (geo%tPeriodic) then
      write(stdOut,*) 'supercell lattice vectors:'
      write(stdOut,*) 'a1:',geo%latVecs(:,1)
      write(stdOut,*) 'a2:',geo%latVecs(:,2)
      write(stdOut,*) 'a3:',geo%latVecs(:,3)
    end if

  end subroutine readGeometry

  !!* Read geometry information for transport calculation
  subroutine readTransportGeometry(root, geom, tp)
    type(fnode), pointer :: root
    type(TGeometry), intent(inout) :: geom
    type(TTransPar), intent(inout) :: tp

    type(fnode), pointer :: pGeom, pDevice, pNode, pTask, pTaskType
    type(string) :: modif
    type(fnode), pointer :: pTmp, field
    type(fnodelist), pointer :: pNodeList
    !type(fnodeList), pointer :: pNodeList
    integer :: ii, contact
    real(dp) :: acc, contactRange(2), sep

    tp%defined = .true.
    tp%tPeriodic1D = .not. geom%tPeriodic
    call getChild(root, "Device", pDevice)
    call getChildValue(pDevice, "AtomRange", tp%idxdevice)
    call getChild(pDevice, "FirstLayerAtoms", pTmp, requested=.false.)
    call readFirstLayerAtoms(pTmp, tp%PL, tp%nPLs, tp%idxdevice)
    if (.not.associated(pTmp)) then
      call setChildValue(pDevice, "FirstLayerAtoms", tp%PL)
    end if
    call getChildren(root, "Contact", pNodeList)
    tp%ncont = getLength(pNodeList)
    if (tp%ncont < 2) then
      call detailedError(pGeom, "At least two contacts must be defined")
    end if
    allocate(tp%contacts(tp%ncont))
    !! Parse contact geometry
    call readContacts(pNodeList, tp%contacts, geom, .true.)

  end subroutine readTransportGeometry

  !> Reads settings for the first layer atoms in principal layers
  subroutine readFirstLayerAtoms(pnode, pls, npl, idxdevice, check)

    type(fnode), pointer, intent(in) :: pnode

    !> Start atoms in the principal layers
    integer, allocatable, intent(out) :: pls(:)

    !> Number of principal layers
    integer, intent(out) :: npl

    !> Atoms range of the device
    integer, intent(in) :: idxdevice(2)

    !> Optional setting to turn on/off check (defaults to on if absent)
    logical, optional, intent(in) :: check


    type(TListInt) :: li
    logical :: checkidx

    checkidx = .true.
    if (present(check)) checkidx = check

    if (associated(pnode)) then
        call init(li)
        call getChildValue(pnode, "", li)
        npl = len(li)
        allocate(pls(npl))
        call asArray(li, pls)
        call destruct(li)
        if (checkidx) then
          if (any(pls < idxdevice(1) .or. &
                  pls > idxdevice(2))) then
             call detailedError(pnode, "First layer atoms must be between " &
               &// i2c(idxdevice(1)) // " and " // i2c(idxdevice(2)) // ".")
          end if
        end if
      else
         npl = 1
         allocate(pls(npl))
         pls = (/ 1 /)
      end if

  end subroutine readFirstLayerAtoms

   !> Read bias information, used in Analysis and Green's function solver
  subroutine readContacts(pNodeList, contacts, geom, upload)
    type(fnodeList), pointer :: pNodeList
    type(ContactInfo), allocatable, dimension(:), intent(inout) :: contacts
    type(TGeometry), intent(in) :: geom
    logical, intent(in) :: upload

    real(dp) :: contactLayerTol, vec(3)
    integer :: ii, jj
    type(fnode), pointer :: field, pNode, pTmp, pWide
    type(string) :: buffer, modif
    type(TListReal) :: fermiBuffer

    do ii = 1, size(contacts)

      contacts(ii)%wideBand = .false.
      contacts(ii)%wideBandDos = 0.0

      call getItem1(pNodeList, ii, pNode)
      call getChildValue(pNode, "Id", buffer, child=pTmp)
      buffer = tolower(trim(unquote(char(buffer))))
      if (len(buffer) > mc) then
        call detailedError(pTmp, "Contact id may not be longer than " &
            &// i2c(mc) // " characters.")
      end if
      contacts(ii)%name = char(buffer)
      if (any(contacts(1:ii-1)%name == contacts(ii)%name)) then
        call detailedError(pTmp, "Contact id '" // trim(contacts(ii)%name) &
            &//  "' already in use")
      end if

      call getChildValue(pNode, "PLShiftTolerance", contactLayerTol, 1e-5_dp, modifier=modif,&
          & child=field)
      call convertByMul(char(modif), lengthUnits, field, contactLayerTol)

      call getChildValue(pNode, "AtomRange", contacts(ii)%idxrange, child=pTmp)
      call getContactVectorII(contacts(ii)%idxrange, geom, ii, pTmp, contactLayerTol, &
                              & contacts(ii)%lattice, contacts(ii)%dir)
      contacts(ii)%length = sqrt(sum(contacts(ii)%lattice**2))

      ! Contact temperatures. Needed
      call getChildValue(pNode, "Temperature", contacts(ii)%kbT,&
                         &0.0_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, contacts(ii)%kbT)

      if (upload) then
        contacts(ii)%potential = 0.d0

        call getChildValue(pNode, "wideBand", contacts(ii)%wideBand, .false.)
        if (contacts(ii)%wideBand) then
          call getChildValue(pNode, 'LevelSpacing', contacts(ii)%wideBandDos, &
                             &0.735_dp, modifier=modif, child=field)
          call convertByMul(char(modif), energyUnits, field,&
                              &contacts(ii)%wideBandDos)
          !WideBandApproximation is defined as energy spacing between levels
          !In the code the inverse value (Density of states) is used
          !Convert the negf input value. Default is 20.e eV
          contacts(ii)%wideBandDos = 1.d0 / contacts(ii)%wideBandDos
        end if
        !call getChildValue(pNode, "FermiLevel", contacts(ii)%eFermi,modifier=modif)
        !call convertByMul(char(modif), energyUnits, pNode, contacts(ii)%eFermi)
        contacts(ii)%eFermi=0.d0
      end if

    enddo

  end subroutine readContacts

      ! Sanity checking of atom ranges and returning contact vector and direction.
  subroutine getContactVectorII(atomrange, geom, id, pContact, plShiftTol, &
      &contactVec, contactDir)
    integer, intent(in) :: atomrange(2)
    type(TGeometry), intent(in) :: geom
    integer, intent(in) :: id
    type(fnode), pointer :: pContact
    real(dp), intent(in) :: plShiftTol
    real(dp), intent(out) :: contactVec(3)
    integer, intent(out) :: contactDir

    integer :: iStart, iStart2, iEnd
    logical :: mask(3)

    !! Sanity check for the atom ranges
    iStart = atomrange(1)
    iEnd = atomrange(2)
    if (iStart < 1 .or. iEnd < 1 .or. iStart > geom%nAtom &
        &.or. iEnd > geom%nAtom .or. iEnd < iStart) then
      call detailedError(pContact, "Invalid atom range '" // i2c(iStart) &
          &// " " // i2c(iEnd) // "', values should be between " // i2c(1) &
          &// " and " // i2c(geom%nAtom) // ".")
    end if
    if (mod(iEnd - iStart + 1, 2) /= 0) then
      call detailedError(pContact, "Nr. of atoms in the contact must be even")
    end if

    ! Determining contact vector
    iStart2 = iStart + (iEnd - iStart + 1) / 2
    contactVec = geom%coords(:,iStart) - geom%coords(:,iStart2)
    if (any(sqrt(sum((geom%coords(:,iStart:iStart2-1) - geom%coords(:,iStart2:iEnd) &
        &- spread(contactVec, dim=2, ncopies=iStart2-iStart))**2, dim=1)) > plShiftTol)) then
      write(stdOut,*) 'coords:', geom%coords(:,iStart)
      write(stdOut,*) 'coords:', geom%coords(:,iStart2)
      write(stdOut,*) 'Contact Vector:', contactVec(1:3)
      write(stdOut,*) iStart,iStart2,iEnd
      write(stdOut,*) 'X:'
      write(stdOut,*) ((geom%coords(1,iStart:iStart2-1) - geom%coords(1,iStart2:iEnd)&
          & - spread(contactVec(1), dim=1, ncopies=iStart2-iStart)))
      write(stdOut,*) 'Y:'
      write(stdOut,*) ((geom%coords(2,iStart:iStart2-1) - geom%coords(2,iStart2:iEnd) &
          & - spread(contactVec(2), dim=1, ncopies=iStart2-iStart)))
      write(stdOut,*) 'Z:'
      write(stdOut,*) ((geom%coords(3,iStart:iStart2-1) - geom%coords(3,iStart2:iEnd) &
          &- spread(contactVec(3), dim=1, ncopies=iStart2-iStart)))
      call error("Contact " // i2c(id) &
          &// " does not consist of two rigidly shifted layers."//new_line('a') &
          &// "Check structure or increase PLShiftTolerance.")
    end if

    contactDir = 0

  end subroutine getContactVectorII

  !> Used to read atomic masses from SK files
  subroutine readSKfiles(child, geo, speciesMass)
    type(fnode), pointer :: child
    type(TGeometry), intent(in) :: geo
    real(dp), dimension(:) :: speciesMass

    type(TOldSKData) :: skData
    type(TListCharLc), allocatable :: skFiles(:)
    type(fnode), pointer :: value, child2
    type(string) :: buffer, buffer2
    character(lc) :: prefix, suffix, separator, elem1, elem2, strTmp, filename
    type(TListString) :: lStr
    integer :: ii, iSp1
    logical :: tLower, tExist

    write(stdOut, "(/, A)") "read atomic masses from sk files..."
    !! Slater-Koster files
    allocate(skFiles(geo%nSpecies))
    do iSp1 = 1, geo%nSpecies
        call init(skFiles(iSp1))
    end do

    call getChildValue(child, "", value)
    call getNodeName(value, buffer)

    select case(char(buffer))
    case ("type2filenames")
      call getChildValue(value, "Prefix", buffer2, "")
      prefix = unquote(char(buffer2))
      call getChildValue(value, "Suffix", buffer2, "")
      suffix = unquote(char(buffer2))
      call getChildValue(value, "Separator", buffer2, "")
      separator = unquote(char(buffer2))
      call getChildValue(value, "LowerCaseTypeName", tLower, .false.)
      do iSp1 = 1, geo%nSpecies
        if (tLower) then
          elem1 = tolower(geo%speciesNames(iSp1))
        else
          elem1 = geo%speciesNames(iSp1)
        end if
        strTmp = trim(prefix) // trim(elem1) // trim(separator) &
            &// trim(elem1) // trim(suffix)
        call append(skFiles(iSp1), strTmp)
        inquire(file=strTmp, exist=tExist)
        if (.not. tExist) then
          call detailedError(value, "SK file with generated name '" &
              &// trim(strTmp) // "' does not exist.")
        end if
      end do
    case default
      call setUnprocessed(value)
      do iSp1 = 1, geo%nSpecies
        strTmp = trim(geo%speciesNames(iSp1)) // "-" &
            &// trim(geo%speciesNames(iSp1))
        call init(lStr)
        call getChildValue(child, trim(strTmp), lStr, child=child2)
        ! We can't handle selected shells here (also not needed I guess)
        if (len(lStr) /= 1) then
          call detailedError(child2, "Incorrect number of Slater-Koster &
              &files")
        end if
        do ii = 1, len(lStr)
          call get(lStr, strTmp, ii)
          inquire(file=strTmp, exist=tExist)
          if (.not. tExist) then
            call detailedError(child2, "SK file '" // trim(strTmp) &
                &// "' does not exist'")
          end if
          call append(skFiles(iSp1), strTmp)
        end do
        call destruct(lStr)
      end do
    end select

    do iSp1 = 1, geo%nSpecies
      call get(skFiles(iSp1), fileName, 1)
      call readFromFile(skData, fileName, .true.)
      deallocate(skData%skHam)
      deallocate(skData%skOver)
      speciesMass(iSp1) = skData%mass
    end do

    do iSp1 = 1, geo%nSpecies
      call destruct(skFiles(iSp1))
    end do
    deallocate(skFiles)

  end subroutine readSKfiles

  subroutine readMasses(value, geo, speciesMass)
    type(fnode), pointer :: value
    type(TGeometry), intent(in) :: geo
    real(dp), dimension(:) :: speciesMass

    type(fnode), pointer :: child, child2
    type(string) :: modif
    integer :: iSp
    character(lc) :: strTmp
    real(dp) :: mass, defmass

    write(stdOut, "(/, A)") "set atomic masses as IUPAC defaults ..."

    do iSp = 1, geo%nSpecies
      defmass = getAtomicMass(trim(geo%speciesNames(iSp)))
      call getChildValue(value, geo%speciesNames(iSp), mass, defmass,&
               &modifier=modif, child= child2)
      speciesMass(iSp) = mass * amu__au
      write(stdOut,*) trim(geo%speciesNames(iSp)),": ", mass, "amu", &
            &SpeciesMass(iSp),"a.u."
    end do

  end subroutine readMasses

  !> K-Points
  subroutine readKPoints(node, geo, tKMesh, tBadIntegratingKPoints)

    !> Relevant node in input tree
    type(fnode), pointer :: node

    !> Geometry structure to be filled
    type(TGeometry), intent(in) :: geo

    !> Whether a KMesh was defined
    logical, intent(out) :: tKMesh

    !> Error check
    logical, intent(out) :: tBadIntegratingKPoints

    type(fnode), pointer :: value1, child
    type(string) :: buffer, modifier
    integer :: ind, ii, jj, kk, n3Tmp(3)
    real(dp) :: coeffsAndShifts(3, 4)
    real(dp) :: r3Tmp(3), r3Tmp2(3), v1(3), v2(3), v3(3)
    type(TListIntR1) :: li1
    type(TListRealR1) :: lr1
    integer, allocatable :: tmpI1(:)
    real(dp), allocatable :: kpts(:,:)
    character(lc) :: errorStr

    ! Assume SCC can has usual default number of steps if needed
    tBadIntegratingKPoints = .false.
    tKMesh = .false.

    ! K-Points
    if (geo%tPeriodic) then
      call getChildValue(node, "KPointsAndWeights", value1, child=child, &
          &modifier=modifier)
      call getNodeName(value1, buffer)

      select case(char(buffer))
      case ("supercellfolding")
        tBadIntegratingKPoints = .false.
        if (len(modifier) > 0) then
          call detailedError(child, "No modifier is allowed, if the &
              &SupercellFolding scheme is used.")
        end if
        call getChildValue(value1, "", coeffsAndShifts)
        if (abs(determinant33(coeffsAndShifts(:,1:3))) - 1.0_dp < -1e-6_dp) then
          call detailedError(value1, "Determinant of the supercell matrix must &
              &be greater than 1")
        end if
        if (any(abs(modulo(coeffsAndShifts(:,1:3) + 0.5_dp, 1.0_dp) - 0.5_dp) &
            &> 1e-6_dp)) then
          call detailedError(value1, "The components of the supercell matrix &
              &must be integers.")
        end if
        call getSuperSampling(coeffsAndShifts(:,1:3), modulo(coeffsAndShifts(:,4), 1.0_dp),&
            & kPoint, kWeight, reduceByInversion=.true.)

        nKPoints = size(kPoint, dim=2)

      case ("kmesh")  
        tBadIntegratingKPoints = .false.
        tKMesh = .true.
        if (len(modifier) > 0) then
          call detailedError(child, "No modifier is allowed, if the &
              &SupercellFolding scheme is used.")
        end if
        call getChildValue(value1, "cells", n3Tmp)
        if (any(n3Tmp<0)) then
           call detailedError(child, "cells must be positive integers")     
        end if
        where (n3Tmp == 0) 
           r3Tmp = 0.0_dp
           r3Tmp2 = real(n3Tmp+1, dp)
        elsewhere
           r3Tmp = 0.5_dp   
           r3Tmp2 = real(n3Tmp, dp)
        end where
        coeffsAndShifts = 0.0_dp
        do ii = 1, 3
          coeffsAndShifts(ii,ii) = r3Tmp2(ii)
        end do
        ! Cell volume calculation: 2D or 3D lattice cases
        v1(:)=geo%latVecs(:,1)  
        v2(:)=geo%latVecs(:,2)  
        v3(:)=geo%latVecs(:,3) 
        if (any(n3Tmp==0)) then
          if (minloc(n3Tmp,1)==1) then    
            cellVolume = norm2(cross3(v2,v3))
          else if (minloc(n3Tmp,1)==2) then
            cellVolume = norm2(cross3(v1,v3))
          else 
            cellVolume = norm2(cross3(v1,v2))
          end if      
          print*,'cell Volume:',cellVolume,'Bohr^2'
        else      
          cellVolume = dot_product(v3,cross3(v1,v2))
          print*,'cell Volume:',cellVolume,'Bohr^3'
        end if 

        call getSuperSampling(coeffsAndShifts(:,1:3), r3Tmp, &
            & kPoint, kWeight, reduceByInversion=.false.)

        call kmesh%initMeshFromFolding(n3Tmp, kPoint)

        r3Tmp(1) = kWeight(1)
        deallocate(kPoint)
        deallocate(kWeight)
        nKPoints = size(kmesh%nodeCoords,2)
        allocate(kPoint(3,nKPoints))
        kPoint = kmesh%nodeCoords
        allocate(kWeight(nKPoints))
        kWeight = r3Tmp(1)/4.0_dp
      
      case ("klines")
        ! probably unable to integrate charge for SCC
        tBadIntegratingKPoints = .true.
        call init(li1)
        call init(lr1)
        call getChildValue(value1, "", 1, li1, 3, lr1)
        if (len(li1) < 1) then
          call detailedError(value1, "At least one line must be specified.")
        end if
        allocate(tmpI1(len(li1)))
        allocate(kpts(3, 0:len(lr1)))
        call asVector(li1, tmpI1)
        call asArray(lr1, kpts(:,1:len(lr1)))
        kpts(:,0) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
        call destruct(li1)
        call destruct(lr1)
        if (any(tmpI1 < 0)) then
          call detailedError(value1, "Interval steps must be greater equal to &
              &zero.")
        end if
        nKPoints = sum(tmpI1)
        if (nKPoints < 1) then
          call detailedError(value1, "Sum of the interval steps must be greater &
              &than zero.")
        end if
        ii = 1
        do while (tmpI1(ii) == 0)
          ii = ii + 1
        end do
        allocate(kPoint(3, nKPoints))
        allocate(kWeight(nKPoints))
        ind = 1
        do jj = ii, size(tmpI1)
          if (tmpI1(jj) == 0) then
            cycle
          end if
          r3Tmp = (kpts(:,jj) - kpts(:,jj-1)) / real(tmpI1(jj), dp)
          do kk = 1, tmpI1(jj)
            kPoint(:,ind) = kpts(:,jj-1) + real(kk, dp) * r3Tmp
            ind = ind + 1
          end do
        end do
        kWeight(:) = 1.0_dp
        if (len(modifier) > 0) then
          select case (tolower(char(modifier)))
          case ("relative")
          case ("absolute")
            kPoint(:,:) =  matmul(transpose(geo%latVecs), kPoint)
            kpts(:,:) = matmul(transpose(geo%latVecs), kpts)
          case default
            call detailedError(child, "Invalid modifier: '" // char(modifier) &
                &// "'")
          end select
        end if
        deallocate(tmpI1)
        deallocate(kpts)

      case (textNodeName)

        ! no idea, but assume user knows what they are doing
        tBadIntegratingKPoints = .false.

        call init(lr1)
        call getChildValue(child, "", 4, lr1, modifier=modifier)
        if (len(lr1) < 1) then
          call detailedError(child, "At least one k-point must be defined.")
        end if
        nKPoints = len(lr1)
        allocate(kpts(4, nKPoints))
        call asArray(lr1, kpts)
        call destruct(lr1)
        if (len(modifier) > 0) then
          select case (tolower(char(modifier)))
          case ("relative")
            continue
          case ("absolute")
            kpts(1:3,:) =  matmul(transpose(geo%latVecs), kpts(1:3,:))
          case default
            call detailedError(child, "Invalid modifier: '" // char(modifier) &
                &// "'")
          end select
        end if
        allocate(kPoint(3, nKPoints))
        allocate(kWeight(nKPoints))
        kPoint(:,:) = kpts(1:3, :)
        kWeight(:) = kpts(4, :)
        deallocate(kpts)
      case default
        call detailedError(value1, "Invalid K-point scheme")
      end select
    end if

  end subroutine readKPoints

  subroutine  readKPointsFile(child)
    type(fnode),  pointer ::  child
    type(string) :: text

    call getFirstTextChild(child, text)
    call readKPointsFile_help(child, char(text))

  end subroutine  readKPointsFile

  subroutine readKPointsFile_help(child,text)
    type(fnode),  pointer ::  child
    character(len=*), intent(in) :: text
    integer :: iStart, iErr=0, ii, iOldStart
    real(dp), dimension(:), allocatable :: tmparray
    real(dp), dimension(:,:), allocatable :: kpts

    iStart = 1
    call getNextToken(text, nKPoints, iStart, iErr)

    allocate(tmparray(4))
    allocate(kpts(4, nKPoints))
    allocate(kPoint(3, nKPoints))
    allocate(kWeight(nKPoints))

    iErr = -2 !TOKEN_ERROR
    iOldStart = iStart
    iStart  = iOldStart

    do ii = 1, nKPoints
      call getNextToken(text, tmparray, iStart, iErr)
      kpts(:, ii) = tmparray(:)
    end do

    do ii = 1, nKPoints
        kPoint(1:3, ii)  = kpts(1:3, ii)
        kWeight(ii)  = kpts(4, ii)
    end do

  end subroutine readKPointsFile_help



  !>  Read DFTB hessian.
  !>
  !> The derivatives matrix must be stored as the following order:
  !>  For the x y z directions of atoms 1..n
  !>    d^2 E        d^2 E       d^2 E       d^2 E        d^2 E
  !>  ---------- + --------- + --------- + ---------- + ---------- +...
  !>  dx_1 dx_1    dy_1 dx_1   dz_1 dx_1   dx_2 dx_1    dy_2 dx_1
  !>
  subroutine readDftbHessian(child, atMass, dynMat)
    type(fnode), pointer :: child
    real(dp), intent(in) :: atMass(:)
    real(dp), allocatable :: dynMat(:,:)

    type(TListRealR1) :: realBuffer
    integer :: iCount, jCount, ii, kk, jj, ll
    integer :: nDerivs

    integer ::  n, j1, j2, fu
    type(fnode), pointer :: child2
    type(string) :: filename
    logical :: texist

    call getChildValue(child, "Filename", filename, "hessian.out")

    inquire(file=trim(char(filename)), exist=texist )
    if (texist) then
      write(stdOut, "(/, A)") "read dftb hessian '"//trim(char(filename))//"'..."
    else
      call detailedError(child,"Hessian file "//trim(char(filename))//" does not exist")
    end if

    nDerivs = 3 * nMovedAtom
    allocate(dynMat(nDerivs,nDerivs))

    open(newunit=fu, file=trim(char(filename)), action='read')
    do ii = 1,  nDerivs
        read(fu,'(4f16.10)') dynMat(1:nDerivs,ii)
    end do

    ! Note: we read the transpose matrix to avoid temporary arrays (ifort warnings).
    ! It should be symmetric or could be symmetrized here
    dynMat = transpose(dynMat)
    ! mass weight the Hessian matrix to get the dynamical matrix
    iCount = 0
    do ii = 1, nMovedAtom
      do kk = 1, 3
        iCount = iCount + 1
        jCount = 0
        do jj = 1, nMovedAtom
          do ll = 1, 3
            jCount = jCount + 1
            dynMat(jCount,iCount) = dynMat(jCount,iCount) &
                & / (sqrt(atMass(ii)) * sqrt(atMass(jj)))
          end do
        end do
      end do
    end do


  close(fu)

  end subroutine readDftbHessian

  subroutine readDynMatrix(child, dynMat)
    type(fnode), pointer :: child
    real(dp), allocatable :: dynMat(:,:)

    type(TListRealR1) :: realBuffer
    integer :: iCount, jCount, ii, kk, jj, ll
    integer :: nDerivs

    nDerivs = 3 * nMovedAtom
    allocate(dynMat(nDerivs,nDerivs))

    !The derivatives matrix must be stored as the following order:

    ! For the x y z directions of atoms 1..n
    !   d^2 E        d^2 E       d^2 E       d^2 E        d^2 E
    ! ---------- + --------- + --------- + ---------- + ---------- +...
    ! dx_1 dx_1    dy_1 dx_1   dz_1 dx_1   dx_2 dx_1    dy_2 dx_1

    write(stdOut, "(/, A)") "read dynamical matrix..."

    call init(realBuffer)
    call getChildValue(child, "", nDerivs, realBuffer)
    if (len(realBuffer)/=nDerivs) then
      call detailedError(child,"wrong number of derivatives supplied:" &
          & // i2c(len(realBuffer)) // " supplied, " &
          & // i2c(nDerivs) // " required.")
    end if
    call asArray(realBuffer, dynMat)
    call destruct(realBuffer)

  end subroutine readDynMatrix

  subroutine readCp2kHessian(child, atMass, DynMat)
    type(fnode), pointer :: child
    real(dp), intent(in) :: atMass(:)
    real(dp), allocatable :: dynMat(:,:)

    type(TListRealR1) :: realBuffer
    integer :: iCount, jCount, ii, kk, jj, ll
    integer :: nDerivs, nBlocks

    real, dimension(:,:), allocatable :: HessCp2k
    integer ::  n, j1, j2,  p,  q, fu
    type(string) :: filename
    logical :: texist

    nDerivs = 3 * nMovedAtom
    allocate(dynMat(nDerivs,nDerivs))

    call getChildValue(child, "Filename", filename, "hessian.cp2k")
    inquire(file=trim(char(filename)), exist=texist )
    if (texist) then
      write(stdOut, "(/, A)") "read cp2k hessian '"//trim(char(filename))//"'..."
    else
      call detailedError(child,"Hessian file "//trim(char(filename))//" does not exist")
    end if

    !The derivatives matrix must be stored as the following order:

    ! For the x y z directions of atoms 1..n
    !   d^2 E        d^2 E       d^2 E       d^2 E        d^2 E
    ! ---------- + --------- + --------- + ---------- + ---------- +...
    ! dx_1 dx_1    dy_1 dx_1   dz_1 dx_1   dx_2 dx_1    dy_2 dx_1

    open(newunit=fu, file=trim(char(filename)), action='read')
    nBlocks = nDerivs/5.0

    allocate(HessCp2k(nDerivs*nBlocks,5))

    do  ii  = 1,  nDerivs*nBlocks
        read(fu,*) HessCp2k(ii,1:5)
    end do

    do ii = 1,  nBlocks
        do  jj  = 1, nDerivs
            p = 1+5*(ii-1)
            q = 5*ii
          dynMat(jj,p:q) =  HessCp2k(jj + nDerivs*(ii-1),1:5)
        end do
    end do

    ! mass weight the Hessian matrix to get the dynamical matrix
    iCount = 0
    do ii = 1, nMovedAtom
      do kk = 1, 3
        iCount = iCount + 1
        jCount = 0
        do jj = 1, nMovedAtom
          do ll = 1, 3
            jCount = jCount + 1
            dynMat(jCount,iCount) = dynMat(jCount,iCount) &
                & / (sqrt(atMass(ii)) * sqrt(atMass(jj)))
          end do
        end do
      end do
    end do

    close(fu)

  end subroutine readCp2kHessian

  ! LAMMPS Matrix already includes masses 
  ! Units: 1/M d^2 E/dR^2 = 1/E d^2 E/dR^2 c^2 -> 1/T^2
  ! 
  ! LAMMPS stores in units of 1/fs^2 
  ! 1/fs^2 = 5.851001385401955e-4 1/au^2
  ! 
  subroutine readLammpsMatrix(child, dynMat)
    type(fnode), pointer :: child
    real(dp), allocatable :: dynMat(:,:)

    type(string) :: filename
    logical :: texist
    integer :: nDerivs, fu, iCount, jCount, jj
    real(dp), parameter :: conv=5.851001385401955e-4 
    real(dp) :: Dyn(3)

    call getChildValue(child, "Filename", filename, "dynmat.dat")

    inquire(file=trim(char(filename)), exist=texist )
    if (texist) then  
      write(stdOut, "(/, A)") "read Lammps Dynamical Matrix '"//trim(char(filename))//"'..."
    else
      call detailedError(child,"file "//trim(char(filename))//" does not exist")
    end if
    
    open(newunit=fu, file=trim(char(filename)), action='read')
    
    nDerivs = 3 * nMovedAtom
    allocate(dynMat(nDerivs,nDerivs))

    ! mass weight the Hessian matrix to get the dynamical matrix
    do iCount = 1, nDerivs
      jCount = 1
      do jj = 1, nMovedAtom
        read(fu,*) Dyn
        dynMat(jCount:jCount+2,iCount) = Dyn * conv
        jCount = jCount + 3
      end do
    end do

    close(fu)

  end subroutine readLammpsMatrix

  ! VASP Matrix in units of eV/Ang**2 
  ! 1 Bohr = 0.529177 Ang; 1 Hartree = 27.211399 eV
  ! 1 eV/Ang = 0.0194468869 a.u. 
  ! 1 eV/Ang**2 = 0.01029084529351100249
  subroutine readVaspHessian(child, atMass, dynMat)
    type(fnode), pointer :: child
    real(dp), intent(in) :: atMass(:)
    real(dp), allocatable :: dynMat(:,:)

    type(string) :: filename
    logical :: texist
    integer :: nDerivs, nAts, iAt, jAt, iCount, jCount, ii, jj, fu
    real(dp), parameter :: conv=0.0102908452_dp 
    real(dp) :: Dyn(3,3)

    call getChildValue(child, "Filename", filename, "dynmat.dat")

    inquire(file=trim(char(filename)), exist=texist )
    if (texist) then  
      write(stdOut, "(/, A)") "read VASP Forces '"//trim(char(filename))//"'..."
    else
      call detailedError(child,"file "//trim(char(filename))//" does not exist")
    end if
    
    open(newunit=fu, file=trim(char(filename)), action='read')

    read(fu,*) nAts, nAts
    if (nAts .ne. nMovedAtom) then
      call detailedError(child, "Number of moving atoms and forces do not agree")
    end if  
    
    nDerivs = 3 * nMovedAtom
    allocate(dynMat(nDerivs,nDerivs))

    ! mass weight the Hessian matrix to get the dynamical matrix
    do ii = 1, nMovedAtom
      do jj = 1, nMovedAtom
        read(fu,*) iAt, jAt 
        iCount = 3*iAt-2
        jCount = 3*jAt-2
        read(fu,*) Dyn
        dynMat(iCount:iCount+2,jCount:jCount+2) = Dyn * conv &
             & / (sqrt(atMass(iAt)) * sqrt(atMass(jAt)))
      end do
    end do
    close(fu)
    ! Impose the acoustic sum rule
    !do ii = 1, nDerivs
    !  dynMat(ii,ii) = -sum(dynMat(ii,1:ii-1))-sum(dynMat(ii,ii+1:nMovedAtom))
    !end do

  end subroutine readVaspHessian

  ! Subroutine removing entries in the Dynamical Matrix.
  ! Not used because identified as a wrong way
  subroutine selectModes()

    integer :: iCount, jCount, ii, jj, kk, ll

    select case ( selTypeModes )
    case(modeEnum%INPLANE)
      iCount = 0
      do ii = 1, nMovedAtom
        do kk = 1, 3
          iCount = iCount + 1
          jCount = 0
          do jj = 1, nMovedAtom
            do ll = 1, 3
              jCount = jCount + 1
              if (mod(iCount,3).eq.0 .or. mod(jCount,3).eq.0) then
                  dynMatrix(jCount,iCount) = 0.0_dp
              end if
            end do
          end do
        end do
      end do
    case(modeEnum%OUTOFPLANE)
      iCount = 0
      do ii = 1, nMovedAtom
        do kk = 1, 3
          iCount = iCount + 1
          jCount = 0
          do jj = 1, nMovedAtom
            do ll = 1, 3
              jCount = jCount + 1
              if (mod(iCount,3).ne.0 .and. mod(jCount,3).ne.0) then
                  dynMatrix(jCount,iCount) = 0.0_dp
              end if
            end do
          end do
        end do
      end do
    end select

  end subroutine selectModes

  !! Reads the Analysis block.
  subroutine readAnalysis(node, geo, pdos, tundos, transpar, atTemperature)
    type(fnode), pointer :: node, pnode
    type(TGeometry), intent(in) :: geo
    type(TPdos), intent(inout) :: pdos
    type(TNEGFTunDos), intent(inout) :: tundos
    type(TTransPar), intent(inout) :: transpar
    real(dp) :: atTemperature, TempRange(2)

    type(fnode), pointer :: val, child, field
    type(string) :: modif
    type(fnodeList), pointer :: children
    logical :: tBadKpoints

    call getChild(node, "TunnelingAndDOS", child, requested=.false.)
    if (associated(child)) then
      if (.not.tTransport) then
        call detailedError(node, "Tunneling requires Transport block")
      end if
      call readTunAndDos(child, geo, tundos, transpar, maxval(transpar%contacts(:)%kbT) )
    endif

    call getChild(node, "Conductance", child, requested=.false.)
    if (associated(child)) then
      if (.not.tTransport) then
        call detailedError(node, "Conductance requires Transport block")
      end if
      call getChildValue(child, "TempRange", TempRange, modifier=modif,&
      & child=field)
      call convertByMul(char(modif), energyUnits, field, TempRange)

      call getChildValue(child, "TempStep", TempStep, modifier=modif,&
      & child=field)
      call convertByMul(char(modif), energyUnits, field, TempStep)

       TempMin = TempRange(1)
       TempMax = TempRange(2)
    endif

  end subroutine readAnalysis

  subroutine readPDOSRegions(children, geo, iAtInregion, regionLabels)
    type(fnodeList), pointer :: children
    type(TGeometry), intent(in) :: geo
    type(TWrappedInt1), allocatable, intent(out) :: iAtInRegion(:)
    character(lc), allocatable, intent(out) :: regionLabels(:)

    integer :: nReg, iReg
    integer, allocatable :: tmpI1(:)
    type(fnode), pointer :: child, child2
    type(string) :: buffer
    character(lc) :: strTmp

    nReg = getLength(children)
    allocate(regionLabels(nReg))
    allocate(iAtInRegion(nReg))
    do iReg = 1, nReg
      call getItem1(children, iReg, child)
      call getChildValue(child, "Atoms", buffer, child=child2, &
          & multiple=.true.)
      call getSelectedAtomIndices(child2, char(buffer), geo%speciesNames, geo%species, tmpI1)
      iAtInRegion(iReg)%data = tmpI1
      write(strTmp, "('region',I0)") iReg
      call getChildValue(child, "Label", buffer, trim(strTmp))
      regionLabels(iReg) = unquote(char(buffer))
    end do

  end subroutine readPDOSRegions

  !!* Read Tunneling and Dos options from analysis block
  !!* tundos is the container to be filled
  !!* ncont is needed for contact option allocation
  subroutine readTunAndDos(root, geo, tundos, transpar, temperature)
    type(fnode), pointer :: root
    type(TGeometry), intent(in) :: geo
    type(TNEGFTunDos), intent(inout) :: tundos
    type(TTransPar), intent(inout) :: transpar
    real(dp), intent(in) :: temperature

    type(fnode), pointer :: pNode, pTmp, field
    type(fnodeList), pointer :: pNodeList
    integer :: ii, jj, ind, ncont, nKT
    real(dp) :: eRange(2), eRangeDefault(2)
    type(string) :: buffer, modif
    type(TWrappedInt1), allocatable :: iAtInRegion(:)
    logical, allocatable :: tDirectionResInRegion(:)
    character(lc), allocatable :: regionLabelPrefixes(:)

    tundos%defined = .true.
    ncont = transpar%ncont
    call getChildValue(root, "Verbosity", tundos%verbose, 51)
    call getChildValue(root, "WriteLDOS", tundos%writeLDOS, .true.)
    call getChildValue(root, "WriteTunn", tundos%writeTunn, .true.)

    ! Default meaningful: eRange= (0..10*kT]
    ! nKT is set to GreensFunction default, i.e. 10
    ! I avoid an explicit nKT option because I find it confusing here
    ! (it makes sense only out of equilibrium)
    ! What matters is w*[nB(w;T1)-nB(w;T2)] that is finite lim w->0
    nKT = 10
    eRangeDefault(1) = 0.0001_dp
    eRangeDefault(2) = nKT * temperature

    call getChildValue(root, "FreqRange", eRange, eRangeDefault, &
         & modifier=modif, child=field)
    call convertByMul(char(modif), energyUnits, field, eRange)
    tundos%emin = eRange(1)
    tundos%emax = eRange(2)

    if (eRange(1).le.0.d0) then
       call detailedError(root, "FreqRange must be > 0")
    end if
    if (eRange(2).lt.eRange(1)) then
       call detailedError(root, "Emax < Emin")
    end if

    call getChildValue(root, "FreqStep", tundos%estep, 1.0e-5_dp,  &
       & modifier=modif, child=field)

    call convertByMul(char(modif), energyUnits, field, tundos%estep)

    ! Terminal currents
    call getChild(root, "TerminalCurrents", pTmp, requested=.false.)

    if (associated(pTmp)) then
      call getChildren(pTmp, "EmitterCollector", pNodeList)
      allocate(tundos%ni(getLength(pNodeList)))
      allocate(tundos%nf(getLength(pNodeList)))
      do ii = 1, getLength(pNodeList)
        call getItem1(pNodeList, ii, pNode)
        call getEmitterCollectorByName(pNode, tundos%ni(ii),&
            & tundos%nf(ii), transpar%contacts(:)%name)
      end do
      call destroyNodeList(pNodeList)
    else
      allocate(tundos%ni(ncont-1) )
      allocate(tundos%nf(ncont-1) )
      call setChild(root, "TerminalCurrents", pTmp)
      ind = 1
      do ii = 1, 1
        do jj = ii + 1, ncont
          call setChildValue(pTmp, "EmitterCollector", &
              &(/ transpar%contacts(ii)%name, transpar%contacts(jj)%name /))
          tundos%ni(ind) = ii
          tundos%nf(ind) = jj
          ind = ind + 1
        end do
      end do
    end if

    call getChild(root, "DeltaModel", pNode)
    call readDeltaModel(pNode, tundos)

    call getChildValue(root, "BroadeningDelta", tundos%broadeningDelta, &
        &0.0_dp, modifier=modif, child=field)
    call convertByMul(char(modif), energyUnits, field, &
        &tundos%broadeningDelta)

    call getChildren(root, "Region", pNodeList)
    call readPDOSRegions(pNodeList, geo, iAtInRegion, regionLabelPrefixes)
    call destroyNodeList(pNodeList)

    call addAtomResolvedRegion(tundos%dosOrbitals, tundos%dosLabels)

    call setTypeOfModes(root, transpar)


    contains

      !! Adds one region with all the orbitals of the atoms in it.
      subroutine addAtomResolvedRegion(iOrbRegion, regionLabels)

        type(TWrappedInt1), allocatable, intent(out) :: iOrbRegion(:)
        character(lc), allocatable, intent(out) :: regionLabels(:)

        integer :: nRegion, nAtomInRegion, iReg, ind, ii, jj, iAt
        integer :: nIndices

        nRegion = size(iAtInRegion)
        allocate(iOrbRegion(nRegion))
        allocate(regionLabels(nRegion))

        do iReg = 1, nRegion
          nAtomInRegion = size(iAtInRegion(iReg)%data)
          nIndices = 3*nAtomInRegion
          allocate(iOrbRegion(iReg)%data(nIndices))
          ind = 1
          do ii = 1, nAtomInRegion
            iAt = iAtInRegion(iReg)%data(ii)
            do jj = 0, 2
              iOrbRegion(iReg)%data(ind) = 3*iAt - 2 + jj
              ind = ind + 1
            end do
          end do
          regionLabels(iReg) = regionLabelPrefixes(iReg)
        end do

      end subroutine addAtomResolvedRegion


  end subroutine readTunAndDos

  ! Get contacts for terminal currents by name
  subroutine getEmitterCollectorByName(pNode, emitter, collector, contactNames)
    type(fnode), pointer :: pNode
    integer, intent(out) :: emitter, collector
    character(len=*), intent(in) :: contactNames(:)

    type(TListString) :: lString
    character(len=mc) :: buffer
    integer :: ind
    logical :: tFound

    call init(lString)
    call getChildValue(pNode, "", lString)
    if (len(lString) /= 2) then
      call detailedError(pNode, "You must provide two contacts")
    end if
    call get(lString, buffer, 1)
    emitter = getContactByName(contactNames, buffer, pNode)
    call get(lString, buffer, 2)
    collector = getContactByName(contactNames, buffer, pNode)
    call destruct(lString)

  end subroutine getEmitterCollectorByName

  ! Getting the contact by name
  function getContactByName(contactNames, contName, pNode) result(contact)
    character(len=*), intent(in) :: contactNames(:)
    character(len=*), intent(in) :: contName
    type(fnode), pointer :: pNode
    integer :: contact

    logical :: tFound

    tFound = .false.
    do contact = 1, size(contactNames)
      tFound = (contactNames(contact) == contName)
      if (tFound) then
        exit
      end if
    end do
    if (.not. tFound) then
      call detailedError(pNode, "Invalid collector contact name '" &
          &// trim(contName) // "'")
    end if

  end function getContactByName

  ! Set model for w-dependent delta in G.F.
  subroutine readDeltaModel(root, tundos)
    type(fnode), pointer :: root
    type(TNEGFTunDos), intent(inout) :: tundos

    type(fnode), pointer :: pValue, pChild, field
    type(string) :: buffer, modif

    call getChildValue(root, "", pValue, child=pChild)
    call getNodeName(pValue, buffer)
    ! Delta is repeated to allow different defaults if needed
    select case (trim(char(buffer)))
    case("deltasquared")
      call getChildValue(pValue, "Delta", tundos%delta, &
          &0.0001_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, tundos%delta)
      tundos%deltaModel=0
    case("deltaomega")
      call getChildValue(pValue, "Delta", tundos%delta, &
          &0.0001_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, tundos%delta)
      tundos%deltaModel=1
    case("mingo")
      ! As in Numerical Heat transfer, Part B, 51:333, 2007, Taylor&Francis.
      ! Here Delta is just a dimensionless scaling factor
      call getChildValue(pValue, "Delta", tundos%delta, 0.0001_dp)
      ! We set a cutoff frequency of 2000 cm^-1.
      call getChildValue(pValue, "Wmax", tundos%wmax, &
          &0.009_dp, modifier=modif, child=field)
      call convertByMul(char(modif), energyUnits, field, tundos%delta)
      tundos%deltaModel=2
      ! If Emax >> Wmax delta becomes negative
      if (tundos%Emax > tundos%Wmax) then
        call detailedError(pValue,"In Mingo model check Wmax <= Emax")
      end if
    case default
      call detailedError(pValue,"Unknown deltaModel "//trim(char(buffer)))
    end select
  end subroutine ReadDeltaModel

  ! Build a simple neighbor list. Currently does not work for periodic systems.
  ! Have to fix this important point
  subroutine buildNeighbourList()

    integer ::  iAtom, jAtom, ii, jj, kk, PL1, PL2
    !* First guess for nr. of neighbors.
    integer, parameter :: nInitNeighbours = 100
    real :: disAtom, dd(3)
    integer :: nAllAtom
    real(dp) :: mCutoff
    real(dp), allocatable :: coords(:,:), cellVec(:,:), rCellVec(:,:)
    integer, allocatable :: iCellVec(:)

    call TNeighbourlist_init(neighbourList, geo%nAtom, nInitNeighbours)

    mCutoff = 1.0_dp * cutoff

    if (geo%tPeriodic) then
      !! Make some guess for the nr. of all interacting atoms
      nAllAtom = int((real(geo%nAtom, dp)**(1.0_dp/3.0_dp) + 3.0_dp)**3)
      call getCellTranslations(cellVec, rCellVec, geo%latVecs, &
            &geo%recVecs2p, mCutoff)
    else
      nAllAtom = geo%nAtom
      allocate(rCellVec(3, 1))
      rCellVec(:, 1) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
    end if

    allocate(coords(3, nAllAtom))
    allocate(img2CentCell(nAllAtom))
    allocate(iCellVec(nAllAtom))

    call updateNeighbourList(coords, img2CentCell, iCellVec, neighbourList, &
          &nAllAtom, geo%coords, mCutoff, rCellVec, .false.)

    deallocate(coords)
    deallocate(iCellVec)
    deallocate(rCellVec)

    allocate(nNeighbour(geo%nAtom))
    nNeighbour(:) = 0

    call getNrOfNeighboursForAll(nNeighbour, neighbourList, mCutoff)

    ! Check PL size with neighbor list
    do iAtom = 1, transpar%idxdevice(2)
      PL1 = getPL(iAtom)
      do jj = 1, nNeighbour(iAtom)
        jAtom = img2CentCell(neighbourList%iNeighbour(jj,iAtom))
        if (jAtom > transpar%idxdevice(2)) cycle
        PL2 = getPL(jAtom)
        if (.not.(PL1.eq.PL2 .or. PL1.eq.PL2+1 .or. PL1.eq.PL2-1)) then
          write(stdOut,*) 'ERROR: PL size inconsistent with cutoff'
          stop
        end if
      end do
    end do

  end subroutine buildNeighbourList
        
  subroutine getMaps(centralCellAtom, nCells, p2s_map, s2p_map)
    integer, intent(in) :: centralCellAtom(:)
    integer, intent(in) :: nCells(3)
    type(TWrappedInt1), allocatable, intent(out) :: p2s_map(:)
    integer, allocatable, intent(out) :: s2p_map(:)
  
    integer :: ip, is, ii, jj, kk, iCentAtom, ind 
    real(dp) :: v1(3), v2(3), v3(3), px(3)    

    print*,'supercell atoms',geo%nAtom
    print*,'primitive cell atoms',size(centralCellAtom)
    print*,centralcellAtom(:)
    v1(:) = geo%latVecs(:,1)/real(nCells(1),dp)  
    v2(:) = geo%latVecs(:,2)/real(nCells(2),dp)  
    v3(:) = geo%latVecs(:,3)/real(nCells(3),dp)
    print*, v1
    print*, v2
    print*, v3
    allocate(s2p_map(geo%nAtom))
    allocate(p2s_map(size(centralCellAtom)))
    do ip = 1, size(centralCellAtom)
      allocate(p2s_map(ip)%data(nCells(1)*nCells(2)*nCells(3)))
    end do
    print*,'build maps'
    do ip = 1, size(centralCellAtom)
      iCentAtom = centralCellAtom(ip)
      ind = 1  
      do ii = -nCells(1), nCells(1)
        do jj = -nCells(2), nCells(2)
          do kk = -nCells(3), nCells(3)
            px(:) = geo%coords(:,iCentAtom)+ii*v1(:)+jj*v2(:)+kk*v3(:)
            ! search the atom in supercell (slow search)
            do is = 1, geo%nAtom
              if (all(abs(geo%coords(:,is) - px(:)) < 1e-2 )) then
                 p2s_map(ip)%data(ind) = is
                 ind = ind + 1
                 s2p_map(is) = ip
                 exit
              end if   
            end do

          end do
        end do  
      end do
    end do

  end subroutine getMaps

  subroutine cutDynMatrix()

    integer :: iAtom, jAtom, jj
    real(dp), allocatable :: dynMat2(:,:)

    allocate(dynMat2(3*nMovedAtom, 3*nMovedAtom))
    dynMat2 = 0.0_dp

    do iAtom = 1, geo%nAtom
       do jj = 1, nNeighbour(iAtom)
          jAtom = img2CentCell(neighbourList%iNeighbour(jj, iAtom))
          if (neighbourList%neighDist2(jj,iAtom) .le. cutoff**2) then
            dynMat2(3*(iAtom-1)+1:3*(iAtom-1)+3, 3*(jAtom-1)+1:3*(jAtom-1)+3) = &
                dynMatrix(3*(iAtom-1)+1:3*(iAtom-1)+3, 3*(jAtom-1)+1:3*(jAtom-1)+3)
            dynMat2(3*(jAtom-1)+1:3*(jAtom-1)+3, 3*(iAtom-1)+1:3*(iAtom-1)+3) = &
                dynMatrix(3*(iAtom-1)+1:3*(iAtom-1)+3, 3*(jAtom-1)+1:3*(jAtom-1)+3)
          end if
       end do
    end do

    dynMatrix = dynMat2

    deallocate(dynMat2)

  end subroutine cutDynMatrix

  function getPL(iAt) result(PL)
    integer, intent(in) :: iAt
    integer :: PL

    integer :: ii

    do ii = 1, transpar%nPLs-1
      if (iAt>=transpar%PL(ii) .and. iAt<transpar%PL(ii+1)) then
        PL = ii
      end if
    end do
    if ( iAt>=transpar%PL(transpar%nPLs) .and. iAt<=transpar%idxdevice(2)) then
      PL = transpar%nPLs
    endif
    if (iAt > transpar%idxdevice(2)) then
      PL = transpar%nPLs + 1
    endif

  end function getPL

  !> select family of modes to analyze and restrict transmission
  subroutine setTypeOfModes(root, transpar)
    type(fnode), pointer :: root
    type(TTransPar), intent(inout) :: transpar

    type(string) :: buffer

    !selecting the type of modes you want to analyze
    call getchildValue(root, "ModeType", buffer, "all")
    select case(trim(char(buffer)))
    case("all")
      selTypeModes = modeEnum%ALLMODES
    case("along-x")
      selTypeModes = modeEnum%XX
    case("along-y")
      selTypeModes = modeEnum%YY
    case("along-z")
      selTypeModes = modeEnum%ZZ
    case("longitudinal")
      selTypeModes = modeEnum%LONGITUDINAL
      call checkTypeOfModes(root, transpar)
    case("transverse")
      selTypeModes = modeEnum%TRANSVERSE
      call checkTypeOfModes(root, transpar)
    case("in-plane")
      selTypeModes = modeEnum%INPLANE
      call checkTypeOfModes(root, transpar)
    case("out-of-plane")
      selTypeModes = modeEnum%OUTOFPLANE
      call checkTypeOfModes(root, transpar)
    case default
      call detailedError(root,"Unknown type of modes")
    end select

    transpar%typeModes = selTypeModes

  end subroutine setTypeOfModes

  !> Check that the geometry orientation is consistent with selTypeModes
  !> Currently only checks that transport direction is along z
  subroutine checkTypeOfModes(root, tp)
    type(fnode), pointer :: root
    type(TTransPar), intent(inout) :: tp

    select case (selTypeModes)
    case(modeEnum%LONGITUDINAL, modeEnum%TRANSVERSE)
      call checkAlongZ(root, tp)
    case(modeEnum%INPLANE, modeEnum%OUTOFPLANE)
      call checkAlongZ(root, tp)
    case default
    end select

  end subroutine checkTypeOfModes

  ! check that transport direction is along z
  subroutine checkAlongZ(root, tp)
    type(fnode), pointer :: root
    type(TTransPar), intent(inout) :: tp

    real(dp) :: contactVec(3)
    logical :: mask(3)
    integer :: ii

    do ii = 1, size(tp%contacts)
      contactVec = tp%contacts(ii)%lattice
      ! Determine to which axis the contact vector is parallel.
      mask = (abs(abs(contactVec) - sqrt(sum(contactVec**2))) < 1.0e-8_dp)
      if (count(mask) /= 1 .or. .not.mask(3)) then
        call detailedError(root,"Transport direction is not along z")
      end if
    end do
  end subroutine checkAlongZ

  subroutine dns2csr(M, Sp, drop)
    real(dp), intent(inout) :: M(:,:)
    type(z_CSR), intent(out) :: Sp
    real(dp), intent(in) :: drop

    integer :: nnz, ierr, ii, jj, nrow, ncol

    nrow = size(M,1); ncol = size(M,2)

    nnz = 0
    do jj = 1, ncol
      do ii = 1, nrow
        if (abs(M(ii,jj)) < drop) then
           M(ii,jj) = 0.0_dp
        else    
           nnz = nnz + 1
        end if   
      end do
    end do
    
    call create(Sp, nrow, ncol, nnz)

    call rzdnscsr(M, nnz, Sp%nzval, Sp%colind, Sp%rowpnt, ierr)
   
    if (ierr /= 0) then
      STOP 'error in dnscsr conversion' 
    end if

  end subroutine dns2csr 
      
  subroutine rzdnscsr(dns, nzmax, a, ja, ia, ierr)
    real(dp), intent(in) :: dns(:,:)
    integer, intent(in) :: nzmax
    complex(dp), intent(out) :: a(:)
    integer, intent(out) :: ia(:), ja(:)
    integer, intent(out) :: ierr

    integer :: nrow, ncol, next, ii, jj
    ierr = 0

    nrow = size(dns,1)
    ncol = size(dns,2)
    next = 1
    ia(1) = 1
    do ii=1,nrow
       do jj=1, ncol 
          if (abs(dns(ii,jj)) .eq. 0.0_dp) cycle 
          if (next .gt. nzmax) then
             ierr = ii
             return
          endif
          ja(next) = jj
          a(next) = dns(ii,jj)
          next = next+1
       end do    
       ia(ii+1) = next
    end do

  end subroutine rzdnscsr    

end module phonons_initphonons
