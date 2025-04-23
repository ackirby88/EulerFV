module c_printer_module
  interface
  subroutine c_printer(mode, array, len, step) bind(C, name="c_printer")
    use iso_c_binding
    implicit none
    
    integer(C_INT) :: mode,step
    real(C_DOUBLE), dimension(*) :: array
    integer(C_INT) :: len
  end subroutine c_printer
  end interface
end module
