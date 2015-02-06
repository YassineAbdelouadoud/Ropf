cdef extern from "powerflow.h" namespace "powerflow":
cdef cppclass Network:
    Network(std::vector<float>, std::vector<float>, std::vector<int>, std::vector<int>)
    int length
    std::vector<std::map<int,vector<int> > > nodes
    std::vector fathers
	std::vector x
	std::vector r
    void Backward(std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> U)
	void Forward(std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> U)
    std::list BackwardForward_Sweep(std::vector<float> U,
	std::vector<float> Pload,std::vector<float> Qload,double prec)

cdef class PyNetwork
    cdef Network *thisptr
    def __cinit__(self,std::vector<float> r, std::vector<float> x, std::vector<int> from, std::vector<int> to):
        self.thisptr = new Network(r, x, from, to)
    def __dealloc__(self):
        del self.thisptr
    def Backward(self, std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> U):
        self.thisptr.Backward(Pload, Qload, P, Q, U)
    def Forward(self, std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> U):
        self.thisptr.Forward(Pload, Qload, P, Q, U)
    def BackwardForward_Sweep(self, std::vector<float> U,
	std::vector<float> Pload,std::vector<float> Qload,double prec):
        self.thisptr.BackwardForward_Sweep(U, Pload, Qload, prec)






    