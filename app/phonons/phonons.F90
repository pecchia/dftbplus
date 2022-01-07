!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

program phonons
  use dftbp_common_accuracy, only : dp, lc
  use dftbp_common_constants, only : Hartree__cm, Bohr__AA, Hartree__J, Hartree__eV
  use dftbp_common_constants, only : hbar, au__fs, pi
  use dftbp_common_environment
  use dftbp_common_globalenv
  use phonons_initphonons
  use dftbp_io_message
  use dftbp_io_taggedoutput
  use dftbp_math_eigensolver, only : heev
  use dftbp_math_simplealgebra, only : invert33
  use dftbp_type_typegeometry
  use phonons_libnegfint
  use phonons_mesh
  use ln_structure
  !use omp_lib
  implicit none

  type(TEnvironment) :: env
  logical :: tInitialized
  real(dp), allocatable :: tunnTot(:,:), ldosTot(:,:), conductance(:,:)
  real(dp), allocatable :: currLead(:)
  real(dp), allocatable :: eigenValues(:,:)
  complex(dp), allocatable :: eigenModes(:,:,:)
  real(dp), allocatable :: displ(:,:,:)
  logical :: twriteLDOS
  logical :: twriteTunn
  type (TTaggedWriter) :: taggedWriter
  integer :: err, nProcs
  type(z_CSR), pointer :: pCsrMat 
  integer :: ii, iMode
  
  pCsrMat => dynMatCsr 

  call initGlobalEnv()
  call printHeader()
  call TEnvironment_init(env)
  call initProgramVariables(env)
  call TTaggedWriter_init(taggedWriter)

  nProcs = 1
  write(stdOut,*)'Computing environment'
#:if WITH_MPI
  !!print*,env%mpi%globalComm%rank, env%mpi%globalComm%size, tIOProc, env%mpi%globalComm%lead
  nProcs = env%mpi%globalComm%size
#:else
  write(stdOut,*)'Not compiled with MPI enabled'
#:endif

  if (.not. tTransport .and. nProcs>1) then
    call error("Mode calculation is not parallel yet. Run just on 1 node")
    call destructProgramVariables()
  end if
        

  if (tCompModes) then
    call ComputeModes(eigenValues, displ)
  end if

  if (tPhonDispersion) then
    call PhononDispersion(eigenValues, eigenModes, taggedWriter)
  end if

  ! print*,'completness check eigenModes'
  ! call completness_check(eigenModes, kmesh, 1)
  ! call completness_check(eigenModes, kmesh, 4)
  ! call completness_check(eigenModes, kmesh, 6)

  if (tDensityOfStates) then
    call densityOfStates(eigenValues, kmesh)
  end if

  if (tRelaxationTimes) then
    call ComputeRelaxationTimes2(eigenValues, eigenModes, kmesh)
  end if      

  if (tTransport) then
    twriteTunn = .true.
    twriteLDOS = .true.
    call negf_init(env, transpar, tundos, tInitialized)
    call init_tun_proj(selTypeModes, geo%nAtom)

    if (.not. tInitialized) then
      call error("libnegf not initialized")
    end if

    call negf_init_str(geo%nAtom, transpar, neighbourList%iNeighbour, nNeighbour, img2CentCell)
   
    call calc_phonon_current(env, pCsrMat, tunnTot, ldosTot, currLead, conductance, &
                        & twriteTunn, twriteLDOS)  

    if (tWriteTagged) then              
      call writeTaggedOut(taggedWriter, tunnTot, ldosTot, conductance)      
    end if

  end if

  call destructProgramVariables()
  call destructGlobalEnv()

  write(stdOut, "(A)") "Program completed successfully"

contains
    
  subroutine printHeader()
    write (stdOut, "(A)") repeat("=", 80)
    write (stdOut, "(30X,A)") "PHONONS" 
    write (stdOut, "(A,/)") repeat("=", 80)
    write (stdOut, "(A)") "" 
    write (stdOut, "(A)") "Version 0.1" 
    write (stdOut, "(A)") "A tool to compute phonon transmission in nanostructures based on Hessians" 
    write (stdOut, "(A)") "Authors: Alessandro Pecchia, Leonardo Medrano Sandonas" 
    write (stdOut, "(A)") "When using this code, please cite this work:" 
    write (stdOut, "(A)") "Leonardo Medrano Sandonas, Rafaei Gutierrez, Alessandro Pecchia,"
    write (stdOut, "(A)") "Alexander Croy, Gianaurelio Cuniberti, Quantum phonon transport in"
    write (stdOut, "(A)") "nanomaterials: combining atomistic with non-equilibrium Green's functions"
    write (stdOut, "(A)") "techniques, Entropy 21, 735 (2019)" 
    write (stdOut, "(A)") "" 
  end subroutine printHeader

  subroutine ComputeModes(eigenValues, displ)
    real(dp), allocatable, intent(out) :: eigenValues(:,:)
    real(dp), allocatable, intent(out) :: displ(:,:,:)

    integer  :: ii, jj, kk, ll, nAtom, iMode, jCount, iAt, iAtMoved, fu
    character(lc) :: lcTmp, lcTmp2
    real(dp), allocatable :: dynMat(:,:)

    nAtom = geo%nAtom
 
    allocate(dynMat(3*nMovedAtom,3*nMovedAtom))
    dynMat = dynMatrix
    allocate(eigenValues(3*nMovedAtom,1))
 
    ! solve the eigenproblem
    if (tPlotModes) then
      write(stdOut,*) 'Computing vibrational frequencies and modes'
      call heev(dynMat,eigenValues(:,1),'U','V')
    else
      write(stdOut,*) 'Computing vibrational frequencies'
      call heev(dynMat,eigenValues(:,1),'U','N')
    end if
 
    ! take square root of modes (allowing for imaginary modes) and print
    eigenValues(:,1) =  sign(sqrt(abs(eigenValues(:,1))),eigenValues(:,1))
    write(stdOut,*)'Vibrational modes (cm-1):'
    do ii = 1, 3 * nMovedAtom
      write(stdOut,'(f8.2)') eigenValues(ii,1)*Hartree__cm
    end do
    write(stdOut,*)
 
    if (tPlotModes) then
      write(stdOut,*)'Plotting eigenmodes:'
      write(stdOut,*) ModesToPlot(:)
      ! scale mode components on each atom by mass and then normalise total mode
      do ii = 1, nModesToPlot
        iMode = ModesToPlot(ii)
        jCount = 0
        do jj = 1, nMovedAtom
          do ll = 1, 3
            jCount = jCount + 1
            dynMat(jCount,iMode) = dynMat(jCount,iMode) &
                & /sqrt(atomicMasses(jj))
          end do
        end do
        dynMat(:,iMode) = dynMat(:,iMode) &
            & / sqrt(sum(dynMat(:,iMode)**2))
      end do
 
     ! Create displacment vectors for every atom in every mode.
      allocate(displ(3*nAtom, nModesToPlot, 1))
      displ(:,:,:) = 0.0_dp
      do iAt = 1, nAtom
        if (any(iMovedAtoms == iAt)) then
          ! Index of atom in the list of moved atoms
          ! iAt = 9  
          ! iMovedAtoms = [1, 4, 8, 9, 10]
          ! abs(iMovedAtoms-iAt) = [8, 5, 1, 0, 1] => iAtMoved = 4 
          iAtMoved = minloc(abs(iMovedAtoms - iAt), 1)
          do ii = 1, nModesToPlot
            iMode = ModesToPlot(ii)
            displ(3*iAt-2:3*iAt, ii, 1) =  dynMat(3*iAtMoved-2:3*iAtMoved, iMode)
          end do
        end if
      end do
 
 
      if (tAnimateModes) then
        do ii = 1, nModesToPlot
          iMode = ModesToPlot(ii)
          write(lcTmp,"('mode_',I0)")iMode
          write(lcTmp2, "(A,A)") trim(lcTmp), ".xyz"
          open(newunit=fu, file=trim(lcTmp2), position="rewind", status="replace")
          do kk = 1, nCycles
            do ll = 1, nSteps
              write(fu,*)nAtom
              write(fu,*)'Eigenmode',iMode,eigenValues(iMode,1)*Hartree__cm,'cm-1'
              do iAt = 1, nAtom
                write(fu,'(A3,T4,3F10.6)') &
                    & geo%speciesNames(geo%species(iAt)), &
                    & (geo%coords(:,iAt)&
                    & + cos(2.0_dp * pi * real(ll) / real(nSteps))&
                    & * displ(3*iAt-2:3*iAt,ii,1)) * Bohr__AA
              end do
            end do
          end do
          close(fu)
        end do
      else
        open(fu, file="modes.xyz", position="rewind", status="replace")
        do ii = 1, nModesToPlot
          iMode = ModesToPlot(ii)
          write(fu,*)nAtom
          write(fu,*)'Eigenmode',iMode,eigenValues(iMode,1)*Hartree__cm,'cm-1'
          if (tXmakeMol) then ! need to account for its non-standard xyz vector
            ! format:
            do iAt = 1, nAtom
              write(fu,'(A3,T4,3F10.6,A,3F10.6)') &
                  & geo%speciesNames(geo%species(iAt)), &
                  & geo%coords(:,iAt)* Bohr__AA, ' atom_vector ',&
                  & displ(2*iAt-2:3*iAt,ii,1)
            end do
          else ! genuine xyz format
            do iAt = 1, nAtom
              write(fu,'(A3,T4,6F10.6)') &
                  & geo%speciesNames(geo%species(iAt)), &
                  & geo%coords(:,iAt)* Bohr__AA, &
                  & displ(:,iAt,ii)
            end do
          end if
        end do
        close(fu)
      end if
      
    end if

    deallocate(dynMat)

  end subroutine ComputeModes

  subroutine PhononDispersion(eigenValues, eigenModes, tWriter)
    real(dp), allocatable :: eigenValues(:,:)
    complex(dp), allocatable :: eigenModes(:,:,:)
    type(TTaggedWriter) :: tWriter

    integer  :: ii, jj, kk, ll, iAtom,  iK, jAtom,  kAtom, iAtMoved
    integer  :: i1, j1, i2, k2, fu, ftag, nrep, iMode
    real(dp) ::  ModKPoint,  ModDeltaR
    character(lc) :: lcTmp, lcTmp2
    complex(dp), dimension(:,:), allocatable :: KdynMatrix
    real(dp), dimension(:), allocatable :: tmpvec
    real(dp) :: latVecs(3,3), invLatt(3,3), norm 
    real(dp) :: DeltaR(3), q(3), qold(3)
    complex(dp), parameter ::    j = (0.d0, 1.d0)
    complex(dp) :: phase
    real(dp) :: unitsConv
      
    call setConversionUnits(unitsConv)

    write(stdOut,*) 'Supercell repetitions:' 
    write(stdOut,*) nCells(1),'x',nCells(2),'x',nCells(3) 

    latVecs(:,1) = geo%latVecs(:,1)/real(nCells(1),dp) 
    latVecs(:,2) = geo%latVecs(:,2)/real(nCells(2),dp) 
    latVecs(:,3) = geo%latVecs(:,3)/real(nCells(3),dp) 

    call invert33(invLatt, latVecs) 
    write(stdOut,*) 'reciprocal lattice vectors:'
    write(stdOut,*) 'b1:',invLatt(1,:)
    write(stdOut,*) 'b2:',invLatt(2,:)
    write(stdOut,*) 'b3:',invLatt(3,:)

    invLatt = invLatt * 2.0_dp * pi

    allocate(KdynMatrix(3*nAtomUnitCell,3*nAtomUnitCell))
    allocate(eigenValues(nModesToPlot, nKPoints))
    if (tPlotModes) then
      allocate(eigenModes(3*nMovedAtom, nModesToPlot, nKPoints))
      eigenModes=(0.0_dp,0.0_dp)
    end if

    write(stdOut,*) 'Computing Phonon Dispersion (units '//trim(outputUnits)//')'
    if (tIOProc) then
      open(newunit=fu, file='phononDispersion.dat', action='write')
      if (tWriteTagged) then
        open(newunit=ftag, file=autotestTag) 
      end if
    end if   

    qold(1) = dot_product(invLatt(:,1), kPoint(:,1))
    qold(2) = dot_product(invLatt(:,2), kPoint(:,1))
    qold(3) = dot_product(invLatt(:,3), kPoint(:,1))
    ModKpoint=0.0_dp

    !do iK = 1, nKPoints     
    !  print*, iK,kPoint(1:3,iK),kweight(iK)
    !end do    
    allocate(tmpvec(3*nAtomUnitCell))

    do iK  = 1, nKPoints
      q(1) = dot_product(invLatt(:,1), kPoint(:,iK))
      q(2) = dot_product(invLatt(:,2), kPoint(:,iK))
      q(3) = dot_product(invLatt(:,3), kPoint(:,iK))

      KdynMatrix = 0.d0
      do  ii = 1, nAtomUnitCell
        iAtom = centralCellAtom(ii)
        do  jj = 1, nAtomUnitCell
          jAtom = centralCellAtom(jj)
          ! This loops over all periodic copies
          do  kk  = 1, size(p2s_map(jj)%data)
            kAtom = p2s_map(jj)%data(kk) 
            DeltaR(:) = geo%Coords(:,kAtom)-geo%Coords(:,jAtom)
            ! Phase factor: exp[i q (Rk - Ri)]
            ! (Rk - Ri) = Rk - Rj + (Rj - Ri) 
            ! since kAtoms = equiv(jAtom) => Rk-Rj is always a lattice translation, R
            ! Rj-Ri = Position of atom j relative to i within the primitive cell 
            ! Note: Bloch states here are defined as:
            !       |a,q> = 1/sqrt(Ncells) Sum_R |a,R> exp[i q (R + Ra)]
            !  <a,q|D|b,q> = 1/Ncells Sum_R,R' <a,0|D|b,R'-R> exp[i q (R'-R)]
            ! With this the dynMat in k-space is: Dab(q) = Sum_R <a,0|H|b,R> exp[i q (Rb-Ra)]
            i1 = 3*(ii-1)+1
            j1 = 3*(jj-1)+1
            i2 = 3*(iAtom-1)+1
            k2 = 3*(kAtom-1)+1
            KdynMatrix(i1:i1+2,j1:j1+2) = KdynMatrix(i1:i1+2,j1:j1+2) &
               & + dynMatrix(i2:i2+2,k2:k2+2)*exp(j*dot_product(q,DeltaR))
          end do
        end do
      end do
      
      ! solve the eigenproblem
      if (tPlotModes) then
        call heev(KdynMatrix,tmpvec,'U','V')

        ! Build supercell eigenmodes with proper phases
        ! normalization factor
        norm = sqrt(real(nCells(1)*nCells(2)*nCells(3)))
        ! actual central cell atom 1 for phase factor   
        do iAtom = 1, geo%nAtom
          if (any(iMovedAtoms == iAtom)) then
            ! map of atom to central cell index for eigenvector    
            i1 = s2p_map(iAtom)
            jAtom = centralCellAtom(i1)
            ! Phase factor w.r.t. central cell atoms: exp[i q (Ri - Rk)] 
            ! Includes the relative cell distance and relative distance in p.c. 
            ! This phase is consistent with the definition of D(q) *see above*
            DeltaR(:) = geo%Coords(:,iAtom)-geo%Coords(:,jAtom)
            phase = exp(j*dot_product(q, DeltaR))
            ! Index of atom in the list of moved atoms
            ! iAtom = 9:  iMovedAtoms = [1, 4, 8, 9, 10]
            ! abs(iMovedAtoms-iAtom) = [8, 5, 1, 0, 1] => iAtMoved = 4 
            iAtMoved = minloc(abs(iMovedAtoms - iAtom), 1)
            do ii = 1, nModesToPlot
              iMode = ModesToPlot(ii)
              eigenModes(3*iAtMoved-2:3*iAtMoved, ii, iK) =  &
                    KdynMatrix(3*i1-2:3*i1, iMode) * phase/norm
            end do
          end if
        end do
      else
        call heev(KdynMatrix,tmpvec,'U','N')
      end if


      ! take square root of modes (allowing for imaginary modes) and print
      tmpvec = sign(sqrt(abs(tmpvec)),tmpvec)
      do ii = 1, nModesToPlot
        iMode = ModesToPlot(ii) 
        eigenValues(ii,iK) =  tmpvec(iMode)
      end do
 
      if (tIOProc) then
        ModKPoint = ModKPoint + sqrt(dot_product(q-qold,q-qold)) 
        qold = q 
        do ii = 1, 3*nAtomUnitCell
          write(fu,*) ModKPoint,  tmpvec(ii)*unitsConv
        end do

        if (tWriteTagged) then
          call tWriter%write(ftag, "kpoint", kPoint(:,iK))
          call tWriter%write(ftag, "bands", tmpvec)
        end if
      end if

    end do
    close(fu)
    deallocate(KdynMatrix)
    deallocate(tmpvec)
        
  end subroutine PhononDispersion

  ! ---------------------------------------------------------------------- 
  subroutine orthonormality_check(Mat)
    complex(dp) :: Mat(:,:)    
        
    integer :: ncol
    real(dp), dimension(:,:), allocatable :: BB

    ncol= size(Mat,2)
    allocate(BB(ncol,ncol))
    BB = abs(matmul(conjg(transpose(Mat(:,:))),Mat(:,:)))
    do ii = 1, ncol
      BB(ii,ii)=abs(BB(ii,ii)-1.0_dp)
    end do 
    if (any(BB > 1e-12)) then
      print*, 'check failed'   
    end if
    deallocate(BB) 
  end subroutine orthonormality_check
  
  
  ! ---------------------------------------------------------------------- 
  subroutine completness_check(eigenModes, kmesh, iaR)
    complex(dp), intent(in) :: eigenModes(:,:,:)
    type(Tmesh), intent(in) :: kmesh
    integer, intent(in) :: iaR

    integer :: iElem, iQ, ibS
    complex(dp) :: cont
    real(dp) :: norm
        
    norm = real(nCells(1)*nCells(2)*nCells(3))*kWeight(1)*4.0_dp

    do ibS = 1, 3*nMovedAtom
      cont = (0.0_dp, 0.0_dp) 
      !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iElem,iQ)
      !$OMP DO REDUCTION(+:cont)
      do iElem = 1, size(kmesh%elem), 4
        iQ = kmesh%elem(iElem)%node(1)%id 
        cont = cont + dot_product(eigenModes(ibS,:,iQ),eigenModes(iaR,:,iQ))*norm
      end do
      !$OMP END DO
      !$OMP END PARALLEL 
      if (ibS == iaR) then
        if (abs(abs(cont)-1.0_dp) > 1e-10) then
           print*,'failed',abs(cont)   
        end if
      else  
        if (abs(cont) > 1e-10) then
           print*,'failed',iaR,ibS,abs(cont)   
        end if
      end if
    end do

  end subroutine completness_check  

  ! ---------------------------------------------------------------------- 
  subroutine densityOfStates(eigenValues, kmesh)
    real(dp), intent(in) :: eigenValues(:,:)
    type(Tmesh), intent(in) :: kmesh

    integer :: nW, ii, fu 
    real(dp), allocatable :: dos(:)
    real(dp) :: Wmin, Wmax, Wq, dW
    real(dp) :: unitsConv

    ! set regular frequency grid of 0.25 THz = 3.8e-5 Hartree 
    dW = 0.0000038_dp
    Wmin = minval(eigenValues)
    Wmax = maxval(eigenValues)
    nW = (Wmax-Wmin)/dW + 1
    allocate(dos(nW))
    call setConversionUnits(unitsConv)
    print*,'Wmin=',Wmin*unitsConv
    print*,'Wmax=',Wmax*unitsConv
    
    ! DENSITY OF STATES --------------------------------------------------
    dos = 0.0_dp
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ii,Wq)
    !$OMP DO
    do ii = 1, nW
      Wq = Wmin + dW * (ii-1) 
      call elem_contributions(kmesh, Wq, dos(ii))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    ! --------------------------------------------------------------------
      
    if (tIOProc) then
      open(newunit=fu, file='DOS.dat', action='write')
      do ii = 1, nW
        write(fu,*) (Wmin+dW*(ii-1))*unitsConv, dos(ii)
      end do
      close(fu)
    end if
    
    deallocate(dos)

  end subroutine densityOfStates

  ! ---------------------------------------------------------------------- 
  ! Relaxation times
  ! tau_ql = pi/(2*w_ql^2) Sum_k,l' |<q,l|V|k,l'>|^2 delta(w_ql-w_kl')
  !
  ! dos = pi Sum_k,l' delta(w_ql-w_kl')
  !
  ! Implementation using linear interpolation on a triangular mesh
  ! 
  subroutine ComputeRelaxationTimes2(eigenValues, eigenModes, kmesh)
    real(dp), intent(in) :: eigenValues(:,:)
    complex(dp), intent(in) :: eigenModes(:,:,:)
    type(Tmesh), intent(in) :: kmesh

    integer :: nW, iQ, iElem, iMode, fu 
    real(dp), allocatable :: tau(:,:)
    real(dp) :: Wq, mytau
    !real(dp) :: latVecs(3,3), invLatt(3,3) 
    complex(dp), allocatable :: V(:,:), Ve(:)
    complex(dp), allocatable :: TM(:,:,:)
    real(dp) :: unitsConv
    integer, external :: omp_get_thread_num
      
    ! Check the mesh is traingular
    do iElem = 1, size(kmesh%elem)
      if (kmesh%elem(iElem)%nsides /= 3) then
         write(stdOut,*) 'ERROR: elements must be triangular'
         stop
      end if    
    end do

    call setConversionUnits(unitsConv)

    ! Get reciprocal lattice vectors 
    !latVecs(:,1) = geo%latVecs(:,1)/real(nCells(1),dp) 
    !latVecs(:,2) = geo%latVecs(:,2)/real(nCells(2),dp) 
    !latVecs(:,3) = geo%latVecs(:,3)/real(nCells(3),dp) 
    !call invert33(invLatt, latVecs) 
    !invLatt = invLatt * 2.0_dp * pi

    ! Perturbation matrix 
    allocate(V(3*nMovedAtom, 3*nMovedAtom))
    V = dynMatrixPert - dynMatrix
    write(stdOut, *) '|V|=',maxval(abs(V))
    write(stdOut, *) '(nK * nModes):', (nKpoints*nModesToPlot)
    
    !write(stdOut, "(A,ES8.4,A)") 'Memory required:', 2*(nKpoints*nModesToPlot/1e4)**2 * 1.6,'GB'  
    allocate(TM(3*nMovedAtom, nModesToPlot, nKpoints))
    write(stdOut, *) 'Compute T-matrix'
    call computeTMatrix(V, kmesh, delta, TM)
    
    allocate(tau(nModesToPlot, nKPoints))
    tau = 0.0_dp

    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iQ,iMode,Wq,mytau), FIRSTPRIVATE(Ve)
    allocate(Ve(3*nMovedAtom))
    !$OMP DO
    ! Loop over |iQ,l>
    do iQ = 1, nKPoints
      !write(stdOut,*) 'Kpoint:',iQ,'/',nKpoints
      do iMode = 1, nModesToPlot
        ! Relaxation times
        Ve = TM(:,iMode,iQ)
        Wq = eigenValues(iMode,iQ)

        ! calculation of element contributions-------------
        call elem_contributions(kmesh, Wq, mytau, Ve)

        tau(iMode,iQ) = mytau/(2.0_dp*Wq**2) 
      end do
    end do
    !$OMP END DO
    deallocate(Ve)
    !$OMP END PARALLEL
    
    
    if (tIOProc) then
      open(newunit=fu, file='RelaxationTimes.dat', action='write')
      ! Relaxations times are written in m^2/s
      ! -> Conversion factor Bohr_AA^2/au__fs * 1e-5
      !    1 AA^2/fs = 1e-20 m^2/1e-15 s = 1e-5 m^2/s
      do iQ = 1, nKPoints
        do iMode = 1, nModesToPlot
          write(fu,*) abs(eigenValues(iMode,iQ))*unitsConv, &
                & tau(iMode,iQ) * cellVolume*Bohr__AA**2/au__fs * 1.0d-5
        end do 
      end do
      close(fu)
    end if   

    deallocate(V)
    deallocate(tau)

  end subroutine ComputeRelaxationTimes2
 
  ! ---------------------------------------------------------------------- 
  ! Contribution of DOS or |<k|V|k'>|^2 over a triangualar mesh
  ! A linear interpolation is used  
  subroutine elem_contributions(kmesh, Wq, scalar, Ve)
    type(TMesh), intent(in) :: kmesh
    real(dp), intent(in) :: Wq
    real(dp), intent(inout) :: scalar
    complex(dp), intent(in), optional :: Ve(:) 

    integer :: jj, iM, iElem, iNode, iK, nElem, P(3)
    real(dp) :: Wk(3), F(3), F1, F2, vol, DD
    complex(dp) :: matel

    nElem = size(kmesh%elem)

    !ModesToPlot, nModesToPlot, eigenValues, eigenModes   are global
    ! calculation of element contributions-------------
    scalar = 0.0_dp
    F = 1.0_dp

    do iM = 1, nModesToPlot
      do iElem = 1, nElem 
        do iNode = 1, 3 
          iK = kmesh%elem(iElem)%node(iNode)%id
          vol = kmesh%elem(iElem)%volume
          Wk(iNode) = eigenValues(iM, iK)
          if (present(Ve)) then
            ! Loop over <iK,l'| -> P = |<iK,l'|V|iQ,l>|^2
            matel = dot_product(eigenModes(:, iM, iK), Ve(:))
            F(iNode) = abs(matel)**2
          end if
        end do
        ! reorder the node labels such that Wk(1) < Wk(2) < Wk(3)
        call sortWk(Wk, P)

        ! Checks the 2 cases:
        !  CASE 1    CASE 2  
        !  2---3     2---3   
        !  |--/      | \/     
        !  | /       | /     
        !  1         1           
        if (Wk(1) < Wq .and. Wq < Wk(2)) then
          DD = (Wq-Wk(1))/((Wk(2)-Wk(1))*(Wk(3)-Wk(1)))   
          F1 = (F(P(2))*(Wq-Wk(1)) + F(P(1))*(Wk(2)-Wq))/(Wk(2)-Wk(1))   
          F2 = (F(P(3))*(Wq-Wk(1)) + F(P(1))*(Wk(3)-Wq))/(Wk(3)-Wk(1))
          scalar = scalar + pi * vol * (F1+F2) * DD
        else if (Wk(2) < Wq .and. Wq < Wk(3)) then
          DD = (Wk(3)-Wq)/((Wk(3)-Wk(1))*(Wk(3)-Wk(2)))   
          F1 = (F(P(3))*(Wq-Wk(2)) + F(P(2))*(Wk(3)-Wq))/(Wk(3)-Wk(2))   
          F2 = (F(P(3))*(Wq-Wk(1)) + F(P(1))*(Wk(3)-Wq))/(Wk(3)-Wk(1))
          scalar = scalar + pi * vol * (F1+F2) * DD
        end if
         
      end do
    end do
  
  end subroutine elem_contributions 

  ! ---------------------------------------------------------------------- 
  ! Sort k-points and P with increasing order of Wk
  subroutine sortWk(Wk, P)
    real(dp), intent(inout) :: Wk(3)
    integer, intent(out) :: P(3)
    
    integer :: ii, il
    P = [1, 2, 3]
    do ii = 1, 3
      do il = ii+1, 3
        if (Wk(ii)>Wk(il)) then
          call swapr(Wk(ii), Wk(il))
          call swapi(P(ii), P(il))
        end if  
      end do
    end do
  end subroutine sortWk
 
  subroutine swapi(a,b)
    integer :: a, b
    integer :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swapi

  subroutine swapr(a,b)
    real(dp) :: a, b
    real(dp) :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swapr

  ! -------------------------------------------------------------------------------
  ! Calculation of the T-matrix
  !
  ! Solve system:  T = V + V g T
  !
  ! T_kk' = V_kk' + Sum_k" V_kk" g(k") T_k"k'
  !
  ! Sum_k" [delta_kk" + Sum_l V_kl delta_lk"/(wq^2 - w_l^2 + i*d*wq) ] T_k"k' = V_kk' 
  !
  subroutine computeTMatrix(V, kmesh, delta, TM)
    complex(dp), intent(in) :: V(:,:)
    type(TMesh), intent(in) :: kmesh
    real(dp), intent(in) :: delta
    complex(dp), intent(inout) :: TM(:,:,:)

    integer :: iQ, iElem, iNode, iMode, iK, iM, iter, niters
    complex(dp), parameter :: j = (0.0_dp, 1.0_dp)
    complex(dp), allocatable :: Ve(:)
    complex(dp) :: matel
    real(dp) :: Wq, Wk, vol, norm

    TM = (0.0_dp,0.0_dp)
    norm = real(nCells(1)*nCells(2)*nCells(3))*kWeight(1)*4.0_dp
    niters = 0  
    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iQ,iMode,iElem,iK,iM,Wq,Wk,iter), &
    !$OMP& FIRSTPRIVATE(Ve)
    allocate(Ve(3*nMovedAtom))
    !$OMP DO
    do iQ = 1, nKPoints
      !print*,iQ,'/',nKPoints 
      do iMode = 1, nModesToPlot
        TM(:,iMode,iQ) = matmul(V, eigenModes(:,iMode,iQ))
        ! T|q,m> = V|q,m> + V g V|q,m> = [I+Vg] V|q,m> 
        ! Vq = V|q,m>
        ! T|q,m> = V|q,m> + Sum_k,m  V|k,m> 1/(wq^2 - wk^2 + i wq d) <k,m|V|q,m>
        ! Completness: Sum_k,m |k,m><k,m| = Id
        Wq = eigenValues(iMode,iQ)
        do iter = 1, niters
          call elem_contributions_ve(kmesh, Wq, Ve, TM(:,iMode,iQ))
          TM(:,iMode,iQ) = TM(:,iMode,iQ) + matmul(V,Ve)
        end do  
      end do
    end do
    !$OMP END DO
    deallocate(Ve)
    !$OMP END PARALLEL

    do iMode = 1, nModesToPlot
      matel = 0.0_dp
      !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iElem,iNode,iQ,vol)
      !$OMP DO REDUCTION(+:matel)
      do iElem = 1, size(kmesh%elem)
        vol = kmesh%elem(iElem)%volume
        do iNode = 1, 3
          iQ = kmesh%elem(iElem)%node(iNode)%id 
          matel = matel + &
              & dot_product(eigenModes(:,iMode,iQ), TM(:,iMode,iQ)) *vol/3.0_dp
        end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      write(stdOut, *) iMode,'  |<ql|V|ql>|=',abs(matel)
    end do

  end subroutine computeTMatrix
  
  ! ---------------------------------------------------------------------- 
  ! Contribution of |k>g(w,k)<k|V|q> over a triangualar mesh
  ! A linear interpolation is used  
  subroutine elem_contributions_ve(kmesh, Wq, Vo, Ve)
    type(TMesh), intent(in) :: kmesh
    real(dp), intent(in) :: Wq
    complex(dp), intent(inout) :: Vo(:)
    complex(dp), intent(in) :: Ve(:) 

    integer :: jj, iM, iElem, iNode, iK, nElem, P(3)
    real(dp) :: Wk(3), vol, DD1, DD2
    complex(dp), allocatable :: F(:,:)
    complex(dp), allocatable :: F12(:), F13(:), F23(:), I(:), R(:)
    complex(dp) :: matel

    nElem = size(kmesh%elem)
    allocate(F(size(Vo),3))
    allocate(F12(size(Vo)))
    allocate(F13(size(Vo)))
    allocate(F23(size(Vo)))
    allocate(I(size(Vo)))
    allocate(R(size(Vo)))

    !ModesToPlot, nModesToPlot, eigenValues, eigenModes   are global
    ! calculation of element contributions-------------
    R = (0.0_dp, 0.0_dp)
    I = (0.0_dp, 0.0_dp)
    Vo = (0.0_dp, 0.0_dp)

    do iM = 1, nModesToPlot
      do iElem = 1, nElem 
        do iNode = 1, 3 
          iK = kmesh%elem(iElem)%node(iNode)%id
          vol = kmesh%elem(iElem)%volume
          Wk(iNode) = eigenValues(iM, iK)
          F(:,iNode) = dot_product(eigenModes(:, iM, iK), Ve(:))*eigenModes(:,iM,iK)
        end do
        ! reorder the node labels such that Wk(1) < Wk(2) < Wk(3)
        call sortWk(Wk, P)

        ! Checks the 2 cases:
        !  CASE 1    CASE 2  
        !  2---3     2---3   
        !  |--/      | \/     
        !  | /       | /     
        !  1         1  
        ! Real and imaginary parts are obtained as in 
        ! Ashraff (1987) J phys C: Sol state phys 20 4823 
        ! This has been generalized here for a complex m.e. F(k)
        ! |K,l'><K,l'|V|Q,l>/(E - Ek + i*delta)
        !
        ! g^r(wq,k) = 1/(wq^2 - wk^2 + 2*i*wq*delta)
        ! g^r(wq>0,k) = 1/(wq^2 - wk^2 + i * delta') = P 1/(wq^2 - wk^2) + i*pi*delta(wq^2-wk^2)
        !   
        !
        ! 2F1 + e1 X1 = F12 + F13
        ! 2F3 + e3 X3 = F23 + F13 
        DD1 = (Wq-Wk(1))/((Wk(2)-Wk(1))*(Wk(3)-Wk(1)))   
        DD2 = (Wk(3)-Wq)/((Wk(3)-Wk(1))*(Wk(3)-Wk(2)))   
        F12(:) = (F(:,P(2))*(Wq-Wk(1)) + F(:,P(1))*(Wk(2)-Wq))/(Wk(2)-Wk(1))   
        F13(:) = (F(:,P(3))*(Wq-Wk(1)) + F(:,P(1))*(Wk(3)-Wq))/(Wk(3)-Wk(1))
        F23(:) = (F(:,P(3))*(Wq-Wk(2)) + F(:,P(2))*(Wk(3)-Wq))/(Wk(3)-Wk(2))   
        if (Wk(1) < Wq .and. Wq < Wk(2)) then
          I(:) = I(:) + pi * vol * (F12(:)+F13(:)) * DD1
        else if (Wk(2) < Wq .and. Wq < Wk(3)) then
          I(:) = I(:) + pi * vol * (F13(:)+F23(:)) * DD2
        end if
        ! The real part contributions
        R(:) = R(:) - (F13(:)+F23(:))*log(abs(Wq-Wk(3))/abs(Wq-Wk(2))) * vol * DD2 &
             &      - (F13(:)+F12(:))*log(abs(Wq-Wk(2))/abs(Wq-Wk(1))) * vol * DD1 &
             &      - 2.0_dp * (F(:,P(3))-F(:,P(1)))/(Wk(3)-Wk(1)) * vol &
             &      - (F(:,P(3))*(Wk(2)-Wk(1))+F(:,P(2))*(Wk(3)-Wk(1))-F(:,P(1))*(Wk(3)-Wk(1))) &
             &      * (Wq-Wk(2))/((Wk(3)-Wk(1))*(Wk(3)-Wk(2))*(Wk(2)-Wk(1))) * vol 
       
        Vo = Vo - cmplx(real(R)-aimag(I), real(I)+aimag(R)) 

      end do
    end do
    deallocate(R,I,F12,F13,F23)
    deallocate(F)
  
  end subroutine elem_contributions_ve

  ! ---------------------------------------------------------------------- 
  subroutine setConversionUnits(unitsConv)
    real(dp), intent(out) :: unitsConv

    select case(trim(outputUnits))  
    case("H")
      unitsConv = 1.0_dp    
    case("eV")
      unitsConv = Hartree__eV
    case("meV")    
      unitsConv = Hartree__eV*100.0_dp
    case("cm")
      unitsConv = Hartree__cm
    case("THz")
      unitsConv = Hartree__J/(hbar*2.0_dp*pi)*1e-12_dp
    case default
      unitsConv = 1.0_dp    
    end select     
  end subroutine setConversionUnits

  ! ---------------------------------------------------------------------- 
  subroutine writeTaggedOut(tWriter, tunnTot, ldosTot, conductance)      
    type(TTaggedWriter) :: tWriter
    real(dp), dimension(:,:) :: tunnTot
    real(dp), dimension(:,:) :: ldosTot
    real(dp), dimension(:,:) :: conductance

    integer :: fu

    open(newunit=fu,file=autotestTag,form="formatted", status="replace")
   
    if (size(tunnTot) > 0) then 
      call tWriter%write(fu,"transmission",tunnTot)
    end if
    if (size(ldosTot) > 0) then 
      call tWriter%write(fu,"PDOS",ldosTot)
    end if
    if (size(conductance) > 0) then 
      call tWriter%write(fu,"conductance",conductance)
    end if

    close(fu)

  end subroutine writeTaggedOut

  ! ---------------------------------------------------------------------------------
  ! Compute the relaxation time using T-Matrix approach
  ! tau_ql(w) = 1/w_ql * Sum_k,l' |<q,l|V|k,l'>|^2 Im[ g^r(w) ] 
  !
  ! Im[g^r(w)] = Im[ 1/(w+i*delta)^2 - w_kl'^2) ]
  !            = 2*w*delta/((w^2-w_kl'^2)^2+4*w^2*delta^2)
  !
  ! lim delta -> 0  Im[g^r]= pi * delta(w^2 - w_kl'^2) 
  ! Take w = w_ql
  !
  ! tau_ql = pi/(2*w_ql^2) Sum_k,l' |<q,l|V|k,l'>|^2 delta(w_ql-w_kl')
  !
  ! Implementation using Lorenzian functions
  !
  subroutine ComputeRelaxationTimes(eigenValues, eigenModes, delta)
    real(dp), intent(in) :: eigenValues(:,:)
    complex(dp), intent(in) :: eigenModes(:,:,:)
    real(dp), intent(in) :: delta

    complex(dp), allocatable :: V(:,:), Ve(:)
    complex(dp) :: matel
    real(dp), allocatable :: tau(:,:), dos(:)
    integer :: ii, jj, iMode, iM, nW, iW, iK, iQ, fu, numthreads
    real(dp) :: delta2, dW2, Wq2, Wmin, Wmax


    Wmin = minval(eigenValues)
    Wmax = maxval(eigenValues)
    nW = (Wmax-Wmin)/delta + 1

    allocate(V(3*nMovedAtom, 3*nMovedAtom))
    allocate(tau(nModesToPlot,nKPoints))
    allocate(dos(nW))
    tau = 0.0_dp

    ! Perturbation matrix 
    V = dynMatrixPert - dynMatrix
    write(stdOut, *) '|V|=',maxval(abs(V))
    do ii = 1, nModesToPlot
      matel = 0.0_dp
      do iQ = 1, nKPoints 
        matel = matel + &
            dot_product(eigenModes(:,ii,iQ), matmul(V,eigenModes(:,ii,iQ)))
      end do
      write(stdOut, *) ii,'. |ee|=',maxval(abs(eigenModes(:,ii,iQ))),'  |<ql|V|ql>|=',abs(matel)
    end do

    delta2 = delta*delta

    write(stdOut, *) 'Main Loop'
    write(stdOut, *) '(nK * nModes)^2:', (nKpoints*nModesToPlot)**2
    !write(stdOut, *) 'OMP parallel:',omp_get_max_threads(),' threads'

    !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(iQ,ii,iMode,iK,jj,iM,iW,matel,dW2,Wq2,Ve)
    allocate(Ve(3*nMovedAtom))
    !$OMP DO
    do iQ = 1, nKPoints
      !write(stdOut,*) 'Kpoint:',iQ,'/',nKpoints
      do iMode = 1, nModesToPlot
        ! DOS for convergence check
        do iW = 1, nW
          dW2 = (Wmin + (iW-1)*delta - eigenValues(iMode,iQ))**2
          dos(iW) = dos(iW) + kWeight(iQ) * delta/(dW2 + delta2)/pi
        end do 
        ! Relaxation times
        Ve = matmul(V, eigenModes(:,iMode,iQ))
        tau(iMode,iQ) = 0.0_dp
        do iK = 1, nKPoints 
          do iM = 1, nModesToPlot
            !write(stdOut,*) 'Kpoint:',iQ,'/',nKpoints,'Mode:',ii,'/',nModesToPlot
            matel = dot_product(eigenModes(:,iM,iK), Ve(:))
            Wq2 = eigenValues(iMode,iQ)**2
            dW2 = (Wq2 - eigenValues(iM,iK)**2)**2
            tau(iMode,iQ) = tau(iMode,iQ) + &
                  & kWeight(iK) * 2.0_dp*delta/(dW2 + 4.0_dp*Wq2*delta2) * abs(matel)**2
          end do 
        end do
      end do
    end do
    !$OMP END DO
    deallocate(Ve)
    !$OMP END PARALLEL

    if (tIOProc) then
      open(newunit=fu, file='RelaxationTimes.dat', action='write')
      do iQ = 1, nKPoints
        do iMode = 1, nModesToPlot
          write(fu,*) abs(eigenValues(iMode,iQ)), tau(iMode,iQ)
        end do 
      end do
      close(fu)
      open(newunit=fu, file='DOS.dat', action='write')
      do iW = 1, nW
        write(fu,*) Wmin+delta*(iW-1), dos(iW)
      end do
      close(fu)
    end if   

    deallocate(V)
    deallocate(tau)
    deallocate(dos)

  end subroutine ComputeRelaxationTimes
  
  ! --------------------------------------------------------------------------
  subroutine ComputeGF(w, delta, Gr)
    real(dp), intent(in) :: w 
    real(dp), intent(in) :: delta   
    complex(dp), allocatable, intent(out) :: Gr(:,:)

    complex(dp), allocatable :: ESH(:,:)
    integer :: ii

    allocate(ESH(3*nMovedAtom, 3*nMovedAtom))
    if (.not.allocated(Gr)) then
      allocate(Gr(3*nMovedAtom, 3*nMovedAtom))
    end if

    ! Assemble [w**2 - H] = [w**2 - H0 - V]
    ESH = -dynMatrixPert
    do ii = 1, 3*nMovedAtom
      ESH(ii,ii) = ESH(ii,ii) + cmplx(w**2, 2.0*w*delta, dp) 
    end do 

    ! Call for the Inversion
    call zinv(Gr, ESH, 3*nMovedAtom) 
  
    deallocate(ESH)

  end subroutine ComputeGF

  ! --------------------------------------------------------------------------
  subroutine zinv(inA,A,n)
    complex(dp), dimension(:,:), intent(out) :: inA
    complex(dp), dimension(:,:), intent(in) :: A
    integer, intent(in) :: n
 
    INTEGER :: ipiv(n), info, i
    complex(dp), dimension(:,:), allocatable  :: LU
 
    allocate(LU(n,n))
 
    LU=A
    call zgetrf( n, n, LU, n, ipiv, info )
 
    if (info.ne.0)  then
       write(*,*)
       write(*,*) 'ERROR in LU factorization (zgetrf)',info
       stop
    end if
 
    inA = (0.0_dp, 0.0_dp)
    do i = 1, n
      inA(i,i) = (1.0_dp, 0.0_dp)
    end do
 
    call zgetrs( 'N', n, n, LU, n, ipiv, inA, n, info )
    if (info.ne.0)  then
       write(*,*)
       write(*,*) 'ERROR in INVERSION (zgetrs)',info
       stop
    end if
 
    deallocate(LU)

  end subroutine zinv


end program phonons
