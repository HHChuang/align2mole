!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                               !
!   Program: Overlap two rigid bodies as much as possible.                                      !
!                                                                                               !
!   Input:                                                                                      !
!       1. structure to be varied; $1, extension: .xyz                                          !
!       2. reference structure; $2, extension: .xyz                                             !
!                                                                                               !
!   Output:                                                                                     !
!       1. Std-out: initial RMSD, final RMSD                                                    !
!       2. modified structure; shifted.$1                                                       !
!                                                                                               !
! History:                                                                                      !
! 2019/11/16, Grace, Thanks for the help form Dr. L.P. Wang, UCD.                               !
!   reference code: https://github.com/leeping/forcebalance/blob/master/src/molecule.py!L692    !
!                                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program main 
    implicit none
    integer(4)  :: NAtoms
    character(len=100),allocatable,dimension(:)  :: coord_atoms
    real(8),allocatable,dimension(:,:)  :: coord_var, coord_ref
    real(8),allocatable,dimension(:,:)  :: coord_var_cen, coord_ref_cen
    real(8) :: rmsd1, rmsd2
    real(8),dimension(3)    :: trans
    real(8),dimension(3,3)  :: R = reshape((/ 0.4694443, -0.27083248, -0.8403998, &
    0.85702105, 0.36878014, 0.35988349, &
    0.21245462, -0.88918557, 0.40523087 /), (/3,3/))
    real(8),allocatable,dimension(:,:)  :: coord_var_shifted

    ! 1. Import two structures from input arguments
        call get_NAtoms(NAtoms)
        allocate(coord_atoms(NAtoms))
        allocate(coord_var(NAtoms,3))
        allocate(coord_ref(NAtoms,3))
        call get_coord_atoms(NAtoms,coord_atoms)
        call get_coord(1,NAtoms,coord_var)
        call get_coord(2,NAtoms,coord_ref)

    ! 2. Calculate the initial RMSD
        call get_rmsd(rmsd1,NAtoms,coord_var,coord_ref)

    ! 3. Move two structures to their centroid; translation
    ! output variables: trans (from coord_ref), coord_var_cen, coord_ref
        allocate(coord_var_cen(NAtoms,3))
        allocate(coord_ref_cen(NAtoms,3))
        call get_centroid(trans,NAtoms,coord_var)
        call trans_coord(coord_var_cen,NAtoms,trans,coord_var)
        call get_centroid(trans,NAtoms,coord_ref)
        call trans_coord(coord_ref_cen,NAtoms,trans,coord_ref)
    
    ! 4. Generate rotation matrix by Kabsch algorithm
        call kabsch(R,NAtoms,coord_var_cen,coord_ref_cen)
        ! write(*,*) R(1,1:3) 

    ! 5. Rotate and translate
        allocate(coord_var_shifted(NAtoms,3))
        call get_shifted(coord_var_shifted,NAtoms,coord_var_cen,R,trans)

    ! 6. Export new structure into file named shifted.$1 (*.xyz)
        call output_struc(NAtoms,coord_atoms,coord_var_shifted)

    ! 7. Std-out initial and final RMSD
        call get_rmsd(rmsd2,NAtoms,coord_var_shifted,coord_ref)
        write(*,'(2(F7.4))') rmsd1,rmsd2

    stop
end program main 

subroutine get_NAtoms(NAtoms)
    implicit none 
    integer(4),intent(out)  :: NAtoms
        character(len=100)  :: struc
        logical :: filestat
    
        call GETARG(1,struc)
        INQUIRE(file=struc,exist=filestat)
        if (filestat) then
            open(10,file=struc,status='old')
        else
            write(*,'(A)') TRIM(struc)//" doesn't exist"
            stop
        end if
        read(10,*) NAtoms
        close(10)
        return
    return 
end subroutine get_NAtoms 

subroutine get_coord_atoms(NAtoms,coord_atoms)
    implicit none
    integer(4),intent(in)   :: NAtoms
    character(len=100),dimension(NAtoms),intent(out) :: coord_atoms
    ! locat variable 
    integer(4)  :: i
    character(len=100)  :: filename, buffer

    call GETARG(1,filename)
    open(10,file=TRIM(ADJUSTL(filename)),status='old',action='read')
    read(10,*) buffer
    read(10,*) buffer
    do i=1,NAtoms
        read(10,*) coord_atoms(i),buffer
    end do
    close(10)
    return
end subroutine get_coord_atoms

subroutine get_coord(order,NAtoms,coord)
    implicit none 
    integer(4), intent(in)  :: order, NAtoms 
    real(8),dimension(NAtoms,3),intent(out) :: coord 
    ! local variable 
    integer(4)  :: i 
    character(len=100)  :: filename,buffer

    call GETARG(order,filename)
    open(10,file=TRIM(ADJUSTL(filename)),status='old',action='read')
    read(10,*) buffer
    read(10,*) buffer
    do i=1,NAtoms
        read(10,*) buffer,coord(i,1:3)
    end do
    close(10)
    return 
end subroutine get_coord

subroutine get_rmsd(rmsd,NAtoms,coord_var,coord_ref)
    implicit none
    real(8),intent(inout)   :: rmsd 
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(in)  :: coord_var,coord_ref
    ! local variable 
    integer(4)  :: i,j
    
    rmsd = 0.0D0 
    ! Fortran is column-major 
    do j = 1, 3
        do i = 1, NAtoms
            rmsd = rmsd + (coord_var(i,j) - coord_ref(i,j))**2
        end do 
    end do
    rmsd = SQRT(rmsd/NAtoms)
    return 
end subroutine get_rmsd 

subroutine get_centroid(trans,NAtoms,coord)
    implicit none 
    integer(4),intent(in)   :: NAtoms 
    real(8),dimension(NAtoms,3),intent(in) :: coord
    real(8),dimension(3),intent(out) :: trans
    ! local variable 
    integer(4)  :: i,j

    trans = 0.0D0 
    do j = 1, 3
        do i = 1, NAtoms 
            trans(j) = trans(j) + coord(i,j)
        end do 
        trans(j) = trans(j) / NAtoms 
    end do
    return 
end subroutine get_centroid

subroutine trans_coord(coord_cen,NAtoms,trans,coord)
    implicit none 
    integer(4),intent(in)   :: NAtoms 
    real(8),dimension(3),intent(in) :: trans
    real(8),dimension(NAtoms,3),intent(in)  :: coord
    real(8),dimension(NAtoms,3),intent(inout)   :: coord_cen 
    ! local variable 
    integer(4)  :: i, j

    coord_cen = 0.0D0 
    do j = 1, 3
        do i = 1, NAtoms 
            coord_cen(i,j) = coord(i,j) - trans(j)
        end do 
    end do 
    return 
end subroutine trans_coord

subroutine get_shifted(coord_shifted,NAtoms,coord,R,trans)
    implicit none 
    integer(4),intent(in)   :: NAtoms 
    real(8),dimension(NAtoms,3),intent(in)  :: coord 
    real(8),dimension(3,3),intent(in)  :: R
    real(8),dimension(3),intent(in)  :: trans
    real(8),dimension(NAtoms,3),intent(out)  :: coord_shifted
    ! local variable 
    integer(4)  :: i, j, k
    
    coord_shifted = 0.0D0 
    ! rotation 
    call get_dot_matrix(coord_shifted,NAtoms,3,3,3,coord,R)

    ! translation 
    do j = 1, 3 
        do i = 1, NAtoms 
            coord_shifted(i,j) = coord_shifted(i,j) + trans(j)
        end do
    end do
    return 
end subroutine get_shifted 

subroutine output_struc(NAtoms,coord_atoms,coord)
    implicit none 
    integer(4),intent(in)   :: NAtoms 
    character(len=100),dimension(NAtoms),intent(in) :: coord_atoms
    real(8),dimension(NAtoms,3),intent(in)  :: coord
    !local variables 
    integer(4)  :: i
    character(len=100)  :: filename, buffer

    call GETARG(1,buffer)
    filename='shifted.'//TRIM(ADJUSTL(buffer))
    write(buffer,*) NAtoms
    open(10,file=TRIM(filename),status='replace')
    write(10,'(A)') TRIM(buffer)
    write(10,'(A)') TRIM(filename)
    do i = 1, NAtoms 
        write(buffer,*) TRIM(coord_atoms(i)),coord(i,1:3)
        write(10,'(2A,3(1X,ES20.10))') TRIM(buffer)
    end do
    close(10)
    return 
end subroutine output_struc

subroutine kabsch(R,NAtoms,coord_var,coord_ref)
    implicit none 
    real(8),dimension(3,3),intent(inout)    :: R
    integer(4),intent(in)   :: NAtoms 
    real(8),dimension(NAtoms,3),intent(in)  :: coord_var,coord_ref
    ! local variables
    integer(4)  :: i,j,k
    real(8),dimension(3,NAtoms) :: coord_var_t
    real(8),dimension(3,3)  :: covar

    ! construct covariance matrix
    call get_matrix_t(NAtoms,coord_var_t,coord_var)
    call get_dot_matrix(covar,3,NAtoms,NAtoms,3,coord_var_t,coord_ref)

    ! SVD TODO:

    ! Transposition of v,wt
    ! right-hand coord
    ! Create Rotation matrix R
    return 
end subroutine kabsch 

subroutine get_dot_matrix(matrix,row1,row2,col1,col2,coord1,coord2)
    implicit none 
    integer(4),intent(in)   :: row1,row2,col1,col2
    real(8),dimension(row1,col1),intent(in) :: coord1
    real(8),dimension(row2,col2),intent(in) :: coord2
    real(8),dimension(row1,col2),intent(inout) :: matrix
    ! local variables
    integer(4)  :: i,j,k

    if ( col1 .ne. row2 ) then 
        write(*,'(A)') 'size problem of matrix, get_dot_matrix()'
        write(*,'(A)') 'stop program'
        stop
    end if

    matrix = 0.0D0 
    do i = 1, row1 
        do j = 1, col2
            do k = 1, col1
                matrix(i,j) = matrix(i,j) + coord1(i,k) * coord2(k,j)
            end do
        end do
    end do 
    return 
end subroutine get_dot_matrix

subroutine get_matrix_t(NAtoms,coord_t,coord)
    implicit none 
    integer(4),intent(in)   :: NAtoms 
    real(8),dimension(NAtoms,3),intent(in)  :: coord 
    real(8),dimension(3,NAtoms),intent(out) :: coord_t
    ! local variables 
    integer(4)  :: i,j 
    do i = 1, NAtoms 
        do j = 1, 3 
            coord_t(j,i) = coord(i,j)
        end do 
    end do 
    return 
end subroutine get_matrix_t