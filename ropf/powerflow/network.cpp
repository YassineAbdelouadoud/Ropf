#include "powerflow.h"

using namespace powerflow;

Network::Network(std::vector<float> r, std::vector<float> x, std::vector<int> from, std::vector<int> to){

length = r.size()
fathers = from
r = r
x = x
nodes = create_node_list(std::vector<int> from, std::vector<int> to)
}

Network::Network(){}

void Network::Backward(std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> ){


	}

void Network::Forward(std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> ){


	}


std::list Network::BackwardForward_Sweep(std::vector<float> U,
	std::vector<float> Pload,std::vector<float> Qload,double prec){

	}

std::vector<std::map<int,vector<int> > > create_node_list(std::vector<int> from, std::vector<int> to){


}
