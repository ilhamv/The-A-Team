
void XML_input
( 
	std::string&                                             simName,
	unsigned long long&                                      nhist,          
	double&                                                  Ecut_off,
	double&                                                  tcut_off,
	unsigned long long&                                      trackTime,
	Source_Bank&                                             Sbank,     
	std::stack  < Particle_t >&                              Pbank,        
	std::vector < std::shared_ptr<Surface_t>   >&            Surface,     
	std::vector < std::shared_ptr<Region_t>    >&            Region,    
	std::vector < std::shared_ptr<Nuclide_t>   >&            Nuclide,   
	std::vector < std::shared_ptr<Material_t>  >&            Material, 
	std::vector < std::shared_ptr<Estimator_t> >&            Estimator,
	std::vector < std::shared_ptr<Distribution_t<double>> >& double_distributions,
  	std::vector < std::shared_ptr<Distribution_t<int>>    >& int_distributions,
  	std::vector < std::shared_ptr<Distribution_t<Point_t>>>& point_distributions
)
{

}

