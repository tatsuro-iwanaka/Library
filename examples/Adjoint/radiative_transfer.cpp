#include<iostream>
#include<aad/App>

#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[])
{
	if(argc != 2)
	{
		std::cout << "Usage: ./radiative_transfer <configuration filename>" << std::endl;

		return 1;
	}

	aad::RadiativeTransfer rtm(argv[1]);

	const auto& geo = rtm.geometry();
	std::vector<aad::core::RadiativeLayer> adj_source(1);
	adj_source[0].resize(geo.Ntheta, geo.M);
	double theta = 0.0;
	theta = theta * std::numbers::pi / 180.0;
	int idx = 0;
	double min = 1.0E100;
	for(int i = 0; i < geo.theta_uh.size(); ++i)
	{
		if(std::abs(theta - geo.theta_uh[i]) < min)
		{
			idx = i;
			min = std::abs(theta - geo.theta_uh[i]);
		}
	}
	std::vector<std::vector<std::vector<double>>> reflectance_top = std::vector<std::vector<std::vector<double>>>(geo.Ntheta, std::vector<std::vector<double>>(geo.Ntheta, std::vector<double>(geo.Nphi, 0.0)));
	reflectance_top[0][idx][0] = 1.0;
	adj_source[0].reflectance_m_top_cos = aad::core::computeFourierSeriesCoefficients(reflectance_top, geo)[0];

	rtm.run_adjoint(adj_source);

	const auto& result = rtm.result_adjoint();

	std::ofstream output("vertical_sensitivity.dat");
	output << "#altitude, dJ/dn, dJ/dq, dJ/d(ln q), dJ/dn, dJ/dq, dJ/d(ln q), ..." << std::endl;
	for(int i = 0; i < result.layers.size(); ++i)
	{
		output << std::setprecision(8) << std::scientific << result.altitude[i];

		for(int j = 0; j < result.layers[i].species.size(); ++j)
		{
			output << std::setprecision(8) << std::scientific << ", " << result.layers[i].species[j].number_density << ", " << result.layers[i].species[j].mixing_ratio << ", " << result.layers[i].species[j].log_mixing_ratio;
		}

		output << std::setprecision(8) << std::scientific << std::endl;
	}

	// std::ofstream output_cloud("vertical_sensitivity_cloud.dat");
	// output_cloud << "#altitude, r_g, sigma_g, r_g, sigma_g, ..." << std::endl;
	// for(int i = 0; i < result.layers.size(); ++i)
	// {
	// 	output_cloud << std::setprecision(8) << std::scientific << result.altitude[i];

	// 	for(int j = 1; j < result.layers[i].species.size(); ++j)
	// 	{
	// 		output_cloud << std::setprecision(8) << std::scientific << ", " << result.layers[i].species[j].lnd_r_g << ", " << result.layers[i].species[j].lnd_sigma_g;
	// 	}

	// 	output_cloud << std::setprecision(8) << std::scientific << std::endl;
	// }

	return 0;
}