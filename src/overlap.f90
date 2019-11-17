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
! 2019/11/16, Grace, Transform python code into Fortran version.                                !
!                                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program main 
    implicit none
    integer(4)  :: NAtoms
    real(8) :: rmsd1, rmsd2
    real(8),dimension(3)    :: trans
    real(8),dimension(3,3)  :: R 
    real(8),allocatable,dimension(:,:)  :: coord_var, coord_ref
    real(8),allocatable,dimension(:,:)  :: coord_var_cen, coord_ref_cen
    real(8),allocatable,dimension(:,:)  :: coord_var_shifted
    character(len=100),allocatable,dimension(:)  :: coord_atoms

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
    write(10,'(A)') TRIM(ADJUSTL(buffer))
    write(10,'(A)') TRIM(ADJUSTL(filename))
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
    real(8),dimension(3,3)  :: covar,v,wt
    real(8),dimension(3)    :: s
    logical   :: d
    real(8) :: det_v, det_wt

    ! construct covariance matrix
    call get_matrix_t(NAtoms,coord_var_t,coord_var)
    call get_dot_matrix(covar,3,NAtoms,NAtoms,3,coord_var_t,coord_ref)

    ! SVD 
    ! https://software.intel.com/sites/products/documentation/doclib/mkl_sa/11/mkl_lapack_examples/dgesvd_ex.f.htm
    call get_svd(v,s,wt,covar)

    ! Transposition of v,wt 
    ! http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html
    call get_det(det_v,3,v)
    call get_det(det_wt,3,wt)
    if ( det_v * det_wt .lt. 0 ) then 
        d = .false.
    else 
        d = .true.
    end if

    ! Right-hand coord
    if (d) then 
        s(3) = - s(3)
        do i = 1, 3
            v(i,3) = - v(i,3)
        end do
    end if

    ! Create Rotation matrix R
    call get_dot_matrix(R,3,3,3,3,v,wt)
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

subroutine get_svd(v,s,wt,covar)
    implicit none 
    real(8),dimension(3,3),intent(in)   :: covar
    real(8),dimension(3),intent(inout)  :: s
    real(8),dimension(3,3),intent(inout)    :: v,wt
    !local variables 
    integer(4)  :: LDA, LDU, LDVT
    integer(4)  :: INFO, LWORK ! local scalars
    integer(4),parameter  :: LWMAX=1000
    real(8),dimension(LWMAX) :: WORK

    LDA = 3
    LDU = 3
    LDVT = 3
    ! query the optimal workspace 
    LWORK = -1 
    call DGESVD( 'All', 'All', 3, 3, covar, LDA, s, v, LDU, wt, LDVT, &
                 WORK, LWORK, INFO )
    LWORK = MIN( LWMAX, INT(WORK(1)))
    ! compute SVD
    call DGESVD( 'All', 'All', 3, 3, covar, LDA, s, v, LDU, wt, LDVT, &
                 WORK, LWORK, INFO )
    ! check for convergence 
    if ( INFO .gt. 0 ) then 
        write(*,*) 'The algorithm computing SVD failed to converge. get_svd()'
        stop
    end if
    return 
end subroutine get_svd

!https://dualm.wordpress.com/2012/01/06/computing-determinant-in-fortran/
subroutine get_det(det,dim,mat)
    implicit none 
    integer(4),intent(in)   :: dim 
    real(8),dimension(dim,dim),intent(in)   :: mat
    real(8),intent(out) :: det 
    ! local variables 
    integer(4)  :: LDA, INFO, i
    integer(4),dimension(dim)   :: IPIV
    real(8) :: sgn
    real(8),dimension(dim,dim)  :: mat_copy

    IPIV = 0
    LDA = dim 
    mat_copy = mat
    call DGETRF(dim,dim,mat_copy,LDA,IPIV,INFO)
    
    det = 1.0D0 
    do i = 1, dim 
        det = det * mat_copy(i,i)
    end do 

    sgn = 1.0D0 
    do i = 1, dim 
        if ( IPIV(i) .ne. i ) then 
            sgn = -sgn 
        end if 
    end do

    det = sgn * det 

    return 
end subroutine get_det