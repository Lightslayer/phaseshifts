!This program requires a compiler adhering to Fortran 2003 standard.
program hartfock_main
  implicit none
  integer :: argc, iarg
  character(len=255) :: argv, input_file, log_file, output_dir, atom_dir
  character(len=255) :: prefix(2)
  interface
    subroutine hartfock(inpfile, logfile, outdir, atdir)
      character(len=255) :: inpfile, logfile, outdir, atdir
    end subroutine
  end interface

  !Set default arguments:
  input_file = 'inputA'
  log_file   = ''
  output_dir = './'
  atom_dir   = '~/atlib/'
  !Accessing file whose name begins with '~/' :
  if(atom_dir(1:2)=='~/')then
    call getenv('HOME',prefix(1))
    if (len_trim(prefix(1)) == 0) then !assume Windows host
        call getenv('HOMEDRIVE',prefix(1))
        call getenv('HOMEPATH',prefix(2))
        prefix(1) = trim(prefix(1))//trim(prefix(2))
    end if
    atom_dir=trim(prefix(1))//atom_dir(2:len_trim(atom_dir))
  endif

  !Check if any arguments are found.
  argc = command_argument_count()
  iarg=1
  !loop across options.
  do while(iarg <= argc)
    call get_command_argument(iarg, argv)
    select case(adjustl(argv))
    case("--help", "-h")
        write(*,*)"Executing hf is stored at ~/bin, working directory is '.' ."
        write(*,*)"==================================================================="
        write(*,*)"--input, -i <inputA> : name of input file; default is inputA ."
        write(*,*)"--log, -l <ilogA> : name of log file; default is ilogA ."
        write(*,*)"--output-dir, -o <output_path> : path to log file; default is '.'"
        write(*,*)"--atom_dir, -a <atom_path> : path to charge density files; default is ~/atlib/"
        write(*,*)"-------------------------------------------------------------------"
        write(*,*)"Please contact Eric Shirley <eric.shirley@nist.gov> for queries"
        write(*,*)"and comments."
        stop
    case("--input", "-i")
      if (iarg+1 <= argc) then
        call get_command_argument(iarg+1, argv)
        input_file = adjustl(argv)
        iarg = iarg + 1
      end if
    case("--log", "-l")
      log_file = 'ilogA'
      if (iarg+1 <= argc) then
        call get_command_argument(iarg+1, argv)
        log_file = adjustl(argv)
        iarg = iarg + 1
      end if
    case("--output_dir", "-o")
      if (iarg+1 <= argc) then
        call get_command_argument(iarg+1, argv)
        output_dir = adjustl(argv)
        iarg = iarg + 1
      end if
    case("--atom_dir", "-a")
      if (iarg+1 <= argc) then
        call get_command_argument(iarg+1, argv)
        atom_dir = adjustl(argv)
        iarg = iarg + 1
      end if
    end select
    iarg = iarg + 1
  end do
  call hartfock(input_file,log_file,output_dir,atom_dir)
  stop
end program hartfock_main
