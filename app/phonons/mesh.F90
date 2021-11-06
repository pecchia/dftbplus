module phonons_mesh
  use dftbp_common_accuracy, only : dp 
  use dftbp_math_simplealgebra, only : cross3
  implicit none

  public :: TMesh, TElem, TNode

  type :: TNode
    integer :: id
    real(dp) :: coords(3)    
  end type TNode
 
  type :: TElem
    integer :: nsides    
    type(TNode), allocatable :: node(:)
    real(dp) :: baricenter(3)
    real(dp) :: volume
    contains
    procedure :: initElem
    procedure :: printElem  
  end type TElem

  type :: TMesh
    type(TElem), allocatable :: elem(:)
    real(dp), allocatable :: nodeCoords(:,:) 
    contains
    procedure :: initMesh      
    procedure :: initMeshFromFolding  
    procedure :: addCoords  
    procedure :: printMesh  
  end type TMesh

  contains

  subroutine initMesh(this, nelements)
    class(TMesh), intent(inout) :: this   
    
    !> number of elements 
    integer, intent(in) :: nelements 

    allocate(this%elem(nelements))

  end subroutine initMesh


  subroutine initElem(this, nnodes, coords, indices)
    class(TElem), intent(inout) :: this   

    !> number of nodes
    integer, intent(in) :: nnodes

    !> coordinates of the nodes  (3, nnodes)
    real(dp), intent(in) :: coords(3,nnodes)
    
    !> indices of the nodes in the array (nnodes)
    integer, intent(in) :: indices(nnodes)

    integer :: ii 
    real(dp) :: v1(3), v2(3), v3(3)

    allocate(this%node(nnodes))
    this%nsides = nnodes

    this%baricenter = 0.0_dp
    do ii = 1, nnodes
       !/print*, 'init node',ii, coords(1:3,ii)
       this%node(ii)%coords(:) = coords(:,ii)
       this%node(ii)%id = indices(ii)
       this%baricenter(1:3) = this%baricenter(1:3) + coords(1:3,ii)
    end do
    this%baricenter = this%baricenter/nnodes

    if (nnodes == 3) then
      ! compute 'volume' as triangle area
      v1(:) = coords(:,2)-coords(:,1)  
      v2(:) = coords(:,3)-coords(:,1)  
      this%volume = norm2(cross3(v1,v2))/2.0_dp
    else if (nnodes == 4) then  
      ! compute volume of tetrahydron
      v1(:) = coords(:,2)-coords(:,1)  
      v2(:) = coords(:,3)-coords(:,1)  
      v3(:) = coords(:,4)-coords(:,1)  
      this%volume = dot_product(v3,cross3(v1,v2))/8.0_dp
    end if  

  end subroutine initElem

  ! This mesh is based on supercellFolding with shift 1/2
  ! 
  !  ncells = 6 x 3  
  !
  ! 1 +---+---+---+---+---+---+
  !   |\ /|\ /|\ /|\ /|\ /|\ /| 
  !   | o | o | o | o | o | o |
  !   |/ \|/ \|/ \|/ \|/ \|/ \| 
  !   +---+---+---+---+---+---+
  !   |\ /|\ /|\ /|\ /|\ /|\ /| 
  !   | o | o | o | o | o | o |
  !   |/ \|/ \|/ \|/ \|/ \|/ \| 
  !   +---+---+---+---+---+---+
  !   |\ /|\ /|\ /|\ /|\ /|\ /| 
  !   | o | o | o | o | o | o |
  !   |/ \|/ \|/ \|/ \|/ \|/ \| 
  ! 0 +---+---+---+---+---+---+
  !   0                       1
  !

  subroutine initMeshFromFolding(this, ncells, centers)
    class(TMesh), intent(inout) :: this   

    !> number of cells along reciprocal vectors
    integer, intent(in) :: ncells(3)
    
    !> cell centers 3x #centers
    real(dp), intent(in) :: centers(:,:)


    integer :: nnodes, nelements, nKPoints, ii, kk, nz(2)
    real(dp) :: coords(3,3), dk(3), volume
    integer :: cnt, indx, indices(3)

    !print*,'init Mesh:',ncells

    where (ncells == 0) 
      dk = 0.0_dp
    elsewhere
      dk = 0.5_dp/real(ncells, dp)
    end where

    !print*,'dk=',dk

    if (any(ncells==0)) then
      ! 2D mesh    
      nelements = size(centers,2)
      call this%initMesh(4*nelements)

      call reduce(ncells, nz) 
      nKPoints = nz(1)*nz(2) + (nz(1)+1)*(nz(2)+1)
      !print*,'nKPoints=',nKPoints
      allocate(this%nodeCoords(3,nKPoints)) 
      this%nodeCoords = 0.0_dp

      cnt = 0
      volume = 0.0_dp
      do ii = 1, nelements
        coords(:,1) = centers(:,ii)
        call this%addCoords(cnt, coords(:,1), indx)
        indices(1) = indx
        coords(:,2) = centers(:,ii) + [ -dk(1), -dk(2), dk(3) ]
        call this%addCoords(cnt, coords(:,2), indx)
        indices(2) = indx
        coords(:,3) = centers(:,ii) + [ dk(1), -dk(2), dk(3) ]
        call this%addCoords(cnt, coords(:,3), indx)
        indices(3) = indx
        !print*, 'init elem',4*ii-3,':', indices
        call this%elem(4*ii-3)%initElem(3,coords,indices)
        volume = volume + this%elem(4*ii-3)%volume

        coords(:,1) = centers(:,ii)
        call this%addCoords(cnt, coords(:,1), indx)
        indices(1) = indx
        coords(:,2) = centers(:,ii) + [ dk(1), -dk(2), dk(3) ]
        call this%addCoords(cnt, coords(:,2), indx)
        indices(2) = indx
        coords(:,3) = centers(:,ii) + [ dk(1),  dk(2), dk(3) ]
        call this%addCoords(cnt, coords(:,3), indx)
        indices(3) = indx
        !print*, 'init elem',4*ii-2,':', indices
        call this%elem(4*ii-2)%initElem(3,coords,indices)
        volume = volume + this%elem(4*ii-2)%volume
 
        coords(:,1) = centers(:,ii)
        call this%addCoords(cnt, coords(:,1), indx)
        indices(1) = indx
        coords(:,2) = centers(:,ii) + [ dk(1),  dk(2), dk(3) ]
        call this%addCoords(cnt, coords(:,2), indx)
        indices(2) = indx
        coords(:,3) = centers(:,ii) + [-dk(1),  dk(2), dk(3) ]
        call this%addCoords(cnt, coords(:,3), indx)
        indices(3) = indx
        !print*, 'init elem',4*ii-1,':', indices
        call this%elem(4*ii-1)%initElem(3,coords,indices)
        volume = volume + this%elem(4*ii-1)%volume
 
        coords(:,1) = centers(:,ii)
        call this%addCoords(cnt, coords(:,1), indx)
        indices(1) = indx
        coords(:,2) = centers(:,ii) + [ -dk(1),  dk(2), dk(3) ]
        call this%addCoords(cnt, coords(:,2), indx)
        indices(2) = indx
        coords(:,3) = centers(:,ii) + [ -dk(1), -dk(2), dk(3) ]
        call this%addCoords(cnt, coords(:,3), indx)
        indices(3) = indx
        !print*, 'init elem',4*ii-0,':', indices
        call this%elem(4*ii-0)%initElem(3,coords,indices)
        volume = volume + this%elem(4*ii-0)%volume
      end do  
      print*,'Total volume (should be 1):', volume
      if (abs(volume-1.0_dp)>1e-5) then
        print*,'Internal algorithmic error:',volume,'/=',1   
        stop     
      end if              
      if (cnt /= nKPoints) then
        print*,'Internal algorithmic error:',cnt,'/=',nKPoints   
        stop     
      end if      
   else
     ! 3D Mesh
     ! 12 tetrahedrons
     nKPoints = ncells(1)*ncells(2)*ncells(3) + (ncells(1)+1)*(ncells(2)+1)*(ncells(3)+1)
   end if

  end subroutine initMeshFromFolding

  
  subroutine reduce(ncells,nzcells)
    integer, intent(in) :: ncells(3)
    integer, intent(out) :: nzcells(2)    
    
    if (ncells(1) == 0) then
      nzcells(1) = ncells(2); nzcells(2) = ncells(3)
    else if (ncells(2) == 0) then
      nzcells(1) = ncells(1); nzcells(2) = ncells(3)
    else if (ncells(3) == 0) then
      nzcells(1) = ncells(1); nzcells(2) = ncells(2)
    end if    

  end subroutine reduce

  subroutine addCoords(this, cnt, newCoords, indx)
    class(TMesh) :: this

    !> counter of the last node added
    integer, intent(inout) :: cnt 

    !> coordinates of the current node
    real(dp), intent(in) :: newCoords(3)
    
    !> index of already existing node or 0
    integer, intent(out) :: indx

    integer :: ii

    ! check if node already exist in the vector and return index
    do ii = 1, cnt 
      if (all(this%nodeCoords(1:3,ii) - newCoords(1:3) == 0.0_dp)) then
        indx = ii
        return
      end if
    end do    
    ! add if vector was not found
    cnt = cnt + 1; 
    indx = cnt 
    this%nodeCoords(1:3,cnt) = newCoords(1:3)

  end subroutine addCoords


  subroutine printMesh(this, outUnit)
    class(TMesh), intent(in) :: this

    !> Output unit
    integer, intent(in) :: outUnit    

    integer :: ii

    write(outUnit,*) 'Mesh. #elements: ',size(this%elem) 
    do ii = 1, size(this%elem)
      write(outUnit,*) 'Element', ii
      call this%elem(ii)%printElem(outUnit)
    end do 

  end subroutine printMesh


  subroutine printElem(this, outUnit)
    class(TElem), intent(in) :: this

    !> Output unit
    integer, intent(in) :: outUnit    

    integer :: ii    
    
    do ii = 1, size(this%node)   
      write(outUnit,*) 'Node.: ', ii
      write(outUnit,*) this%node(ii)%coords(:)
    end do  
  end subroutine printElem


end module phonons_mesh 
