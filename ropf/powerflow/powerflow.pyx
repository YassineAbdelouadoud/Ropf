# distutils: language = c++
# distutils: sources = ropf/powerflow/network.cpp
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.list cimport list


cdef extern from "powerflow.h" namespace "powerflow":
    cdef cppclass Network:
        Network(vector[float], vector[float], vector[int], vector[int]) except +
        int length
        vector[map[int,vector[int]]] nodes
        vector fathers
        vector x
        vector r
        void Backward(vector[float] Pload,vector[float] Qload,
        vector[float] P,vector[float] Q,vector[float] U)
        void Forward(vector[float] Pload,vector[float] Qload,
        vector[float] P,vector[float] Q,vector[float] U)
        list BackwardForward_Sweep(vector[float] U,
        vector[float] Pload,vector[float] Qload,double prec)

cdef class PyNetwork:
    cdef Network *thisptr
    def __cinit__(self,vector[float] r, vector[float] x, vector[int] fro, vector[int] to):
        self.thisptr = new Network(r, x, fro, to)
    def __dealloc__(self):
        del self.thisptr
    def Backward(self, vector[float] Pload,vector[float] Qload,
	vector[float] P,vector[float] Q,vector[float] U):
        self.thisptr.Backward(Pload, Qload, P, Q, U)
    def Forward(self, vector[float] Pload,vector[float] Qload,
	vector[float] P,vector[float] Q,vector[float] U):
        self.thisptr.Forward(Pload, Qload, P, Q, U)
    def BackwardForward_Sweep(self, vector[float] U,
	vector[float] Pload,vector[float] Qload,double prec):
        self.thisptr.BackwardForward_Sweep(U, Pload, Qload, prec)






    