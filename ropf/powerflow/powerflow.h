#include <vector>
#include <map>
#include <list>

namespace powerflow {
    class Network {
    int length;

	public :
	std::vector<std::map<int, std::vector<int> > > nodes;
	std::vector<int> fathers;
	std::vector<float> x;
	std::vector<float> r;

	Network(std::vector<float> r, std::vector<float> x, std::vector<int> fro, std::vector<int> to);
	~Network();
	void Backward(std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> U);
	void Forward(std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> U);
    std::list< std::vector<float> > BackwardForward_Sweep(std::vector<float> U,
	std::vector<float> Pload,std::vector<float> Qload, double prec);
 	};
}