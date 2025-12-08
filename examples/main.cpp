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
            layer_omegas[i] = 0.9; // 下層 (Bottom)
        }
		else
		{
            layer_omegas[i] = 0.1; // 上層 (Top)
        }
	}
	
	std::vector<aad::core::RadiativeLayer> layers(n_layer);
	auto geo = aad::core::generateGeometryGaussRadau(n_theta, n_theta * 4 - 1, (n_theta * 4 - 2 / 2));

	for(int i = 0; i < layers.size(); ++i)
	{
		layers[i] = initializeIsotropicLayer(layer_taus[i], layer_omegas[i], geo);
	}

	auto result = aad::core::computeAtmosphere(layers, geo);
	std::cout << result.reflectance_m_top_cos[0](0, 0) << ", " << result.optical_thickness << std::endl;

	// テスト用：有限差分法で計算
	for(int i = 0; i < n_layer; ++i)
	{
		std::vector<aad::core::RadiativeLayer> layers_p(n_layer);
		std::vector<aad::core::RadiativeLayer> layers_m(n_layer);

		double dtau = layers[i].optical_thickness * 1.0E-3;

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

    // 「ナディア反射率に対する感度」を知りたいので、その要素の勾配を 1.0 にする
    aad::core::RadiativeLayer adj_result = result; // 構造をコピー
    
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

	// 5. Chain Rule (初期層の勾配 -> 全層のTauへの勾配)
    std::cout << "\nSensitivity (dJ / dTau_total):" << std::endl;
    for(int i = 0; i < n_layer; ++i)
    {
        double n_doubling_pow = std::pow(2.0, layers[i].n_doubling);
        double omega = layer_omegas[i];

        // A. 倍加計算の内部で tau が直接使われている部分からの勾配 (exp(-tau/mu)など)
        double grad_tau_init = adj_layers[i].optical_thickness;

        // B. 初期化時の行列 (R, T) からの勾配
        // initializeIsotropicLayer で R = C * tau としているので、
        // dJ/dtau = dJ/dR * dR/dtau = dJ/dR * (R / tau)
        // ※ R, T は tau に比例する
        double grad_from_matrices = 0.0;
        
        for(int u=0; u<n_theta; ++u)
		{
            for(int v=0; v<n_theta; ++v)
			{
                // dR / dtau_init = omega / (4 * mu * mu')
                double deriv_coeff = omega / (4.0 * geo.mu(u) * geo.mu(v));
                
                // 各行列成分からの寄与を足し合わせる
                grad_from_matrices += adj_layers[i].reflectance_m_top_cos[0](u, v) * deriv_coeff;
                grad_from_matrices += adj_layers[i].reflectance_m_bottom_cos[0](u, v) * deriv_coeff;
                grad_from_matrices += adj_layers[i].transmittance_m_top_cos[0](u, v) * deriv_coeff;
                grad_from_matrices += adj_layers[i].transmittance_m_bottom_cos[0](u, v) * deriv_coeff;
            }
        }

        // tau_init に対する全勾配
        double total_grad_init = grad_tau_init + grad_from_matrices;
        
        // tau_total に対する勾配に変換
        // tau_init = tau_total / 2^N  =>  d/d_total = d/d_init * (1/2^N)
        double grad_total = total_grad_init / n_doubling_pow;

        std::cout << "Layer " << i << " (Top=" << (n_layer-1) << "): " << std::scientific << std::setprecision(8) << grad_total << std::endl;
    }

	return 0;
}