module m_trajectory
 
use m_precision

implicit none

type trajectory_t
    integer :: frame_size = 0
    integer :: header_size = 0
    integer :: num_frames = 0
    character(len=:), allocatable :: fn
    integer :: file_id = 0
    logical :: isopen = .false.

    contains
        procedure :: create       => traj_create
        procedure :: open         => traj_open
        procedure :: clear        => traj_clear
        procedure :: close        => traj_close
        procedure :: read         => traj_read
        procedure :: append_frame => traj_append_frame
        procedure :: write_frame  => traj_write_frame
        generic   :: init         => create, open
end type trajectory_t
 
contains

!******************************************************************************

subroutine traj_create(this, fn, na)
    !!  Creates a `trajectory_t` object with a new underlying file named `fn`.  If
    !!  `fn` already exists, it will be truncated.  The file `fn` is opened for both
    !!  reading and writing.

    class(trajectory_t), intent(out) :: this
    character(len=*), intent(in) :: fn
    integer, intent(in) :: na
    integer :: file_id
    integer :: i

    !Frame components: nblks, coordinates
    !It is assumed that there is only one molecule and coordinates are w.r.t.
    !the molecule c.o.m. 
    this%frame_size = sizeof_int + 3*na*sizeof_real

    ! Representation of header:
    !     * 1 int for `header_size`
    !     * 1 int for `frame_size`

    this%header_size = 2*sizeof_int

    !Create trajectory file
    open(newunit=file_id, file=fn, access='stream', form='unformatted', &
            action='readwrite', status='replace')

    this%fn = fn
    this%file_id = file_id
    this%num_frames = 0
    this%isopen = .true.

    write(this%file_id) this%header_size
    write(this%file_id) this%frame_size

    end subroutine

!******************************************************************************

subroutine traj_open(this, fn, mode)
    !!  Creates a `trajectory_t` object with a prexisting underlying file named
    !!  `fn`.  If `fn` does not exist, an error will be generated. If
    !!  `mode` == 'rw', the file `fn` is opened for both reading and writing.
    !!  If `mode` == 'r', the file `fn` is opened only for reading.
    !""
    class(trajectory_t), intent(out) :: this
    character(len=*),    intent(in) :: fn
    character(len=*),    intent(in) :: mode
    integer(ip_long) :: file_size
    integer :: i
    
    this%fn = fn

    !Readwrite mode
    if (mode == 'rw') then
        open(newunit=this%file_id, file=this%fn, access='stream', &
            form='unformatted', action='readwrite', status='old')
    !Readonly mode
    else if (mode == 'r') then
        open(newunit=this%file_id, file=this%fn, access='stream', &
            form='unformatted', action='read', status='old')
    end if

    !Get size of the file
    inquire(unit=this%file_id, size=file_size)

    read(this%file_id) this%header_size
    read(this%file_id) this%frame_size

    !Integer division to find `num_frames`
    this%num_frames = int((file_size-this%header_size)/this%frame_size, ip) 
    this%isopen = .true.

    end subroutine

!******************************************************************************

subroutine traj_clear(this)
   !! After a call to this subroutine, all memory within `this` is deallocated,
   !! all components of `this` are reset to zero, and the underlying file is
   !! closed (if open).

    class(trajectory_t), intent(in out) :: this

    call this%close()

    if (allocated(this%fn)) deallocate(this%fn)

    this%file_id = 0
    this%header_size = 0
    this%frame_size = 0
    this%num_frames = 0

    end subroutine

!******************************************************************************

subroutine traj_close(this)
    !! Closes the underlying file of a `trajectory_t`.

    class(trajectory_t), intent(in out) :: this

    if (this%isopen) then
        close(this%file_id)
        this%isopen = .false.
    end if

    end subroutine

!******************************************************************************

subroutine traj_read(this, ver, iframe, istage, nblks, nts, coordinates, ierr)
    !! Read from an open trajectory.

    class(trajectory_t), intent(in) :: this
    integer, intent(in)  :: ver
        !!Trajectory file version. 0: brushpak-0.5, 1: mcpak-0.1, 2: mcpak-0.2
    integer, intent(in)  :: iframe
    integer, intent(out) :: istage
    integer, intent(out) :: nblks
    integer, intent(out) :: nts
    integer, intent(out) :: ierr
    real(rp), dimension(:,:), intent(out) :: coordinates
    real(rp), dimension(3) :: com_
    integer(ip_long) :: offset
    integer :: i

    ierr = 0

    if (iframe > this%num_frames) then
        write(*,*) '`iframe`', iframe, 'out of bounds for `num_frames`', this%num_frames
        ierr = 1
        return
    end if

    offset = this%frame_size
    offset = (iframe-1)*offset + this%header_size + 1

    if (ver==0) then
        nblks = 0
        read(this%file_id, pos=offset) istage, nts, com_, coordinates
    else if (ver==1) then
        istage = 0
        read(this%file_id, pos=offset) nblks, nts, coordinates
    else if (ver==2) then
        istage = 0; nts = 0
        read(this%file_id, pos=offset) nblks, coordinates
    end if

    end subroutine

!******************************************************************************

subroutine traj_append_frame(this, nblks, coordinates)
    !! Write a frame to an open trajectory

    class(trajectory_t), intent(in out) :: this
    integer, intent(in) :: nblks
    real(rp), dimension(:,:), intent(in) :: coordinates
    integer :: iframe

    iframe = this%num_frames + 1

    call traj_write_frame(this, iframe, nblks, coordinates)

    end subroutine

!******************************************************************************

subroutine traj_write_frame(this, iframe, nblks, coordinates)
    !! Write a frame to an open trajectory

    class(trajectory_t), intent(in out) :: this
    integer, intent(in) :: iframe
    integer, intent(in) :: nblks
    real(rp), dimension(:,:), intent(in) :: coordinates
    integer(ip_long) :: offset

    ! Added one to handle appending a frame
    if (iframe > this%num_frames+1) then
        write(*,*) '`iframe`', iframe, 'out of bounds for `num_frames`', &
            (this%num_frames+1)
        return
    end if

    offset = this%frame_size
    offset = (iframe-1)*offset + this%header_size + 1

    write(this%file_id, pos=offset) nblks, coordinates

    if (iframe == this%num_frames + 1) this%num_frames = this%num_frames + 1

    end subroutine

!******************************************************************************

end module m_trajectory
