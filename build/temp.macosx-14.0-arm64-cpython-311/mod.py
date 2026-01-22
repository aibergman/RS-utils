from __future__ import print_function, absolute_import, division
import _mod
import f90wrap.runtime
import logging
import numpy

class Rsutil_Lib(f90wrap.runtime.FortranModule):
    """
    Module rsutil_lib
    
    
    Defined at /Users/andersb/Jobb/RS_tools/source/rsutil_lib.f90 lines 1-76
    
    """
    @staticmethod
    def process_lattice_impl():
        """
        process_lattice_impl()
        
        
        Defined at /Users/andersb/Jobb/RS_tools/source/rsutil_lib.f90 lines 19-76
        
        
        -----------------------------------------------------------------
         Read the lattice namelist from lattice.nml using lattice_io.
        -----------------------------------------------------------------
        """
        _mod.f90wrap_rsutil_lib__process_lattice_impl()
    
    @property
    def new_basis(self):
        """
        Element new_basis ftype=real(rp) pytype=float
        
        
        Defined at /Users/andersb/Jobb/RS_tools/source/rsutil_lib.f90 line 14
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _mod.f90wrap_rsutil_lib__array__new_basis(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            new_basis = self._arrays[array_handle]
        else:
            new_basis = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _mod.f90wrap_rsutil_lib__array__new_basis)
            self._arrays[array_handle] = new_basis
        return new_basis
    
    @new_basis.setter
    def new_basis(self, new_basis):
        self.new_basis[...] = new_basis
    
    @property
    def new_label(self):
        """
        Element new_label ftype=character(len=8) pytype=str
        
        
        Defined at /Users/andersb/Jobb/RS_tools/source/rsutil_lib.f90 line 15
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _mod.f90wrap_rsutil_lib__array__new_label(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            new_label = self._arrays[array_handle]
        else:
            new_label = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _mod.f90wrap_rsutil_lib__array__new_label)
            self._arrays[array_handle] = new_label
        return new_label
    
    @new_label.setter
    def new_label(self, new_label):
        self.new_label[...] = new_label
    
    @property
    def new_atomic_number(self):
        """
        Element new_atomic_number ftype=integer pytype=int
        
        
        Defined at /Users/andersb/Jobb/RS_tools/source/rsutil_lib.f90 line 16
        
        """
        array_ndim, array_type, array_shape, array_handle = \
            _mod.f90wrap_rsutil_lib__array__new_atomic_number(f90wrap.runtime.empty_handle)
        if array_handle in self._arrays:
            new_atomic_number = self._arrays[array_handle]
        else:
            new_atomic_number = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                    f90wrap.runtime.empty_handle,
                                    _mod.f90wrap_rsutil_lib__array__new_atomic_number)
            self._arrays[array_handle] = new_atomic_number
        return new_atomic_number
    
    @new_atomic_number.setter
    def new_atomic_number(self, new_atomic_number):
        self.new_atomic_number[...] = new_atomic_number
    
    @property
    def i(self):
        """
        Element i ftype=integer  pytype=int
        
        
        Defined at /Users/andersb/Jobb/RS_tools/source/rsutil_lib.f90 line 17
        
        """
        return _mod.f90wrap_rsutil_lib__get__i()
    
    @i.setter
    def i(self, i):
        _mod.f90wrap_rsutil_lib__set__i(i)
    
    def __str__(self):
        ret = ['<rsutil_lib>{\n']
        ret.append('    new_basis : ')
        ret.append(repr(self.new_basis))
        ret.append(',\n    new_label : ')
        ret.append(repr(self.new_label))
        ret.append(',\n    new_atomic_number : ')
        ret.append(repr(self.new_atomic_number))
        ret.append(',\n    i : ')
        ret.append(repr(self.i))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

rsutil_lib = Rsutil_Lib()

