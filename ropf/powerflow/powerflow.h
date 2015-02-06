namespace powerflow {
    class Network {
    int length;

	public :
	std::vector<std::map<int,vector<int> > > nodes;
	std::vector fathers;
	std::vector x;
	std::vector r;

	Network(std::vector<float> r, std::vector<float> x, std::vector<int> from, std::vector<int> to);
	~Network();
	void Backward(std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> U);
	void Forward(std::vector<float> Pload,std::vector<float> Qload,
	std::vector<float> P,std::vector<float> Q,std::vector<float> U);
    std::list BackwardForward_Sweep(std::vector<float> U,
	std::vector<float> Pload,std::vector<float> Qload,double prec);

 }