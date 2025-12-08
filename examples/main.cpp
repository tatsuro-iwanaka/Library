#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <aad/Core>

aad::core::RadiativeLayer initializeIsotropicLayer(double target_tau, double omega, const aad::core::Geometry& geo, double initial_tau = 1.0E-8)
{
	aad::core::RadiativeLayer layer;
	layer.optical_thickness = target_tau;
	layer.enable_atmospheric_emission = false;
	layer.is_surface = false;
	layer.n_doubling = 0;

	while(layer.optical_thickness > initial_tau)
	{
		layer.optical_thickness *= 0.5;
		layer.n_doubling++;
	}

	double tau = layer.optical_thickness;

	layer.source_up = Eigen::VectorXd::Zero(geo.Ntheta);
	layer.source_down = Eigen::VectorXd::Zero(geo.Ntheta);
	layer.reflectance_m_top_cos.resize(geo.M + 1);
	layer.reflectance_m_bottom_cos.resize(geo.M + 1);
	layer.transmittance_m_top_cos.resize(geo.M + 1);
	layer.transmittance_m_bottom_cos.resize(geo.M + 1);
	
	layer.reflectance_m_top_sin.resize(geo.M + 1);
	layer.reflectance_m_bottom_sin.resize(geo.M + 1);
	layer.transmittance_m_top_sin.resize(geo.M + 1);
	layer.transmittance_m_bottom_sin.resize(geo.M + 1);

	for(int m=0; m <= geo.M; ++m)
	{
		layer.reflectance_m_top_cos[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.reflectance_m_bottom_cos[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.transmittance_m_top_cos[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.transmittance_m_bottom_cos[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);

		layer.reflectance_m_top_sin[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.reflectance_m_bottom_sin[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.transmittance_m_top_sin[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.transmittance_m_bottom_sin[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
	}

	for (int i = 0; i < geo.Ntheta; ++i)
	{
		for (int j = 0; j < geo.Ntheta; ++j)
		{
			double mu_i = geo.mu(i);
			double mu_j = geo.mu(j);
			
			// Reflection
			double r_val = (omega * tau) / (4.0 * (mu_i * mu_j));
			layer.reflectance_m_top_cos[0](i, j) = r_val;
			layer.reflectance_m_bottom_cos[0](i, j) = r_val;

			// Transmission
			double t_val = (omega * tau) / (4.0 * (mu_i * mu_j));
			layer.transmittance_m_top_cos[0](i, j) = t_val;
			layer.transmittance_m_bottom_cos[0](i, j) = t_val;
		}
	}

	return layer;
}

int main(void)
{
	int n_layer = 20;
	int n_theta = 10;
	std::vector<double> layer_taus(n_layer, 0.2);
	std::vector<double> layer_omegas(n_layer);

	for(int i = 0; i < n_layer; ++i)
	{
		if (i < n_layer / 2)
		{
			layer_omegas[i] = 0.9;
		}
		else
		{
			layer_omegas[i] = 0.1;
		}
	}
	
	std::vector<aad::core::RadiativeLayer> layers(n_layer);
	auto geo = aad::core::generateGeometryGaussRadau(n_theta, n_theta * 4 - 1, (n_theta * 4 - 2) / 2);

	for(int i = 0; i < layers.size(); ++i)
	{
		layers[i] = initializeIsotropicLayer(layer_taus[i], layer_omegas[i], geo);
	}

	auto result = aad::core::computeAtmosphere(layers, geo);

	std::cout << "Central Finite Difference" << std::endl;
	for(int i = 0; i < n_layer; ++i)
	{
		std::vector<aad::core::RadiativeLayer> layers_p(n_layer);
		std::vector<aad::core::RadiativeLayer> layers_m(n_layer);

		double dtau = layers[i].optical_thickness * 1.0E-4;

		for(int j = 0; j < layers.size(); ++j)
		{
			if(i == j)
			{
				layers_p[j] = initializeIsotropicLayer(layer_taus[j] + dtau, layer_omegas[j], geo);
				layers_m[j] = initializeIsotropicLayer(layer_taus[j] - dtau, layer_omegas[j], geo);
			}
			else
			{
				layers_p[j] = initializeIsotropicLayer(layer_taus[j], layer_omegas[j], geo);
				layers_m[j] = initializeIsotropicLayer(layer_taus[j], layer_omegas[j], geo);
			}
		}

		auto result_p = aad::core::computeAtmosphere(layers_p, geo);
		auto result_m = aad::core::computeAtmosphere(layers_m, geo);

		double dj = (result_p.reflectance_m_top_cos[0](0, 0) - result_m.reflectance_m_top_cos[0](0, 0)) / (2.0 * dtau);

		std::cout << "Layer " << i << " (Top=" << (n_layer-1) << "): " << std::scientific << std::setprecision(8) << dj << std::endl;
	}

	aad::core::RadiativeLayer adj_result = result;
	
	adj_result.optical_thickness = 0.0;
	adj_result.source_up.setZero();
	adj_result.source_down.setZero();
	for(int m=0; m<=geo.M; ++m)
	{
		adj_result.reflectance_m_top_cos[m].setZero();
		adj_result.reflectance_m_bottom_cos[m].setZero();
		adj_result.transmittance_m_top_cos[m].setZero();
		adj_result.transmittance_m_bottom_cos[m].setZero();
		adj_result.reflectance_m_top_sin[m].setZero();
		adj_result.reflectance_m_bottom_sin[m].setZero();
		adj_result.transmittance_m_top_sin[m].setZero();
		adj_result.transmittance_m_bottom_sin[m].setZero();
	}

	adj_result.reflectance_m_top_cos[0](0, 0) = 1.0;
	auto adj_layers = aad::core::computeAtmosphere_adjoint(layers, geo, adj_result);

	std::cout << "Adjoint" << std::endl;
	for(int i = 0; i < n_layer; ++i)
	{
		double n_doubling_pow = std::pow(2.0, layers[i].n_doubling);
		double omega = layer_omegas[i];

		double grad_tau_init = adj_layers[i].optical_thickness;

		double grad_from_matrices = 0.0;
		
		for(int u=0; u<n_theta; ++u)
		{
			for(int v=0; v<n_theta; ++v)
			{
				double deriv_coeff = omega / (4.0 * geo.mu(u) * geo.mu(v));
				
				grad_from_matrices += adj_layers[i].reflectance_m_top_cos[0](u, v) * deriv_coeff;
				grad_from_matrices += adj_layers[i].reflectance_m_bottom_cos[0](u, v) * deriv_coeff;
				grad_from_matrices += adj_layers[i].transmittance_m_top_cos[0](u, v) * deriv_coeff;
				grad_from_matrices += adj_layers[i].transmittance_m_bottom_cos[0](u, v) * deriv_coeff;
			}
		}

		double total_grad_init = grad_tau_init + grad_from_matrices;
		
		double grad_total = total_grad_init / n_doubling_pow;

		std::cout << "Layer " << i << " (Top=" << (n_layer-1) << "): " << std::scientific << std::setprecision(8) << grad_total << std::endl;
	}

	return 0;
}