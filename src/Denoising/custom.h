#ifndef CUSTOM_H
#define CUSTOM_H

#include<string>
#include<vector>

using namespace std;


typedef double num_t;

//Uniform Laplacian
const string uniform_key_word = "uniform";

const num_t uniform_num_iterations_default = 10;
const num_t uniform_num_iterations_min = 1;
const num_t uniform_num_iterations_max = 100;
const num_t uniform_num_iterations_delta = 1;

const string uniform_parameter_order = "<number of iterations>";

const vector<num_t>
    uniform_parameters_default = {uniform_num_iterations_default};
const vector<num_t>
    uniform_parameters_min = {uniform_num_iterations_min};
const vector<num_t>
    uniform_parameters_max = {uniform_num_iterations_max};
const vector<num_t>
    uniform_parameters_delta = {uniform_num_iterations_delta};

//Taubin
const string taubin_key_word = "taubin";

const num_t taubin_num_iterations_default = 40;
const num_t taubin_num_iterations_min = 1;
const num_t taubin_num_iterations_max = 100;
const num_t taubin_num_iterations_delta = 1;

const num_t taubin_lambda_default = 0.5;
const num_t taubin_lambda_min = 0.4;
const num_t taubin_lambda_max = 0.6;
const num_t taubin_lambda_delta = 0.01;

const num_t taubin_mu_default = -0.53;
const num_t taubin_mu_min = -0.63;
const num_t taubin_mu_max = -0.43;
const num_t taubin_mu_delta = 0.01;

const string taubin_parameter_order = "<number of iterations> <lambda> <mu>";

const vector<num_t>
    taubin_parameters_default = {taubin_num_iterations_default,taubin_lambda_default,taubin_mu_default};
const vector<num_t>
    taubin_parameters_min = {taubin_num_iterations_min,taubin_lambda_min,taubin_mu_min};
const vector<num_t>
    taubin_parameters_max = {taubin_num_iterations_max,taubin_lambda_max,taubin_mu_max};
const vector<num_t>
    taubin_parameters_delta = {taubin_num_iterations_delta,taubin_lambda_delta,taubin_mu_delta};

//HC Laplacian
const string hc_key_word = "hc";

const num_t hc_num_iterations_default = 40;
const num_t hc_num_iterations_min = 1;
const num_t hc_num_iterations_max = 100;
const num_t hc_num_iterations_delta = 1;

const num_t hc_alpha_default = 0.0;
const num_t hc_alpha_min = 0.0;
const num_t hc_alpha_max = 1.0;
const num_t hc_alpha_delta = 0.01;

const num_t hc_beta_default = 0.5;
const num_t hc_beta_min = 0.0;
const num_t hc_beta_max = 1.0;
const num_t hc_beta_delta = 0.01;

const string hc_parameter_order = "<number of iterations> <alpha> <beta>";

const vector<num_t>
    hc_parameters_default = {hc_num_iterations_default,hc_alpha_default,hc_beta_default};
const vector<num_t>
    hc_parameters_min = {hc_num_iterations_min,hc_alpha_min,hc_beta_min};
const vector<num_t>
    hc_parameters_max = {hc_num_iterations_max,hc_alpha_max,hc_beta_max};
const vector<num_t>
    hc_parameters_delta = {hc_num_iterations_delta,hc_alpha_delta,hc_beta_delta};

//Bilateral
const string bilateral_key_word = "bilateral";

const num_t bilateral_num_iterations_default = 50;
const num_t bilateral_num_iterations_min = 1;
const num_t bilateral_num_iterations_max = 100;
const num_t bilateral_num_iterations_delta = 1;

const string bilateral_parameter_order = "<number of iterations>";

const vector<num_t>
    bilateral_parameters_default = {bilateral_num_iterations_default};
const vector<num_t>
    bilateral_parameters_min = {bilateral_num_iterations_min};
const vector<num_t>
    bilateral_parameters_max = {bilateral_num_iterations_max};
const vector<num_t>
    bilateral_parameters_delta = {bilateral_num_iterations_delta};

//Fast and effective
const string fast_key_word = "fast";

const num_t fast_num_normal_iterations_default = 20;
const num_t fast_num_normal_iterations_min = 1;
const num_t fast_num_normal_iterations_max = 100;
const num_t fast_num_normal_iterations_delta = 1;

const num_t fast_num_vertex_iterations_default = 10;
const num_t fast_num_vertex_iterations_min = 1;
const num_t fast_num_vertex_iterations_max = 100;
const num_t fast_num_vertex_iterations_delta = 1;

const num_t fast_threshold_T_default = 0.5;
const num_t fast_threshold_T_min = 0.1;
const num_t fast_threshold_T_max = 0.9;
const num_t fast_threshold_T_delta = 0.01;

const string fast_parameter_order = "<number of normal iterations> <number of vertex iterations> <threshold_T>";

const vector<num_t>
    fast_parameters_default = {fast_num_normal_iterations_default,fast_num_vertex_iterations_default,fast_threshold_T_default};
const vector<num_t>
    fast_parameters_min = {fast_num_normal_iterations_min,fast_num_vertex_iterations_min,fast_threshold_T_min};
const vector<num_t>
    fast_parameters_max = {fast_num_normal_iterations_max,fast_num_vertex_iterations_max,fast_threshold_T_max};
const vector<num_t>
    fast_parameters_delta = {fast_num_normal_iterations_delta,fast_num_vertex_iterations_delta,fast_threshold_T_delta};


//Bilateral Normal
const string bilateralnormal_key_word = "bilateralnormal";

const num_t bilateralnormal_num_normal_iterations_default = 20;
const num_t bilateralnormal_num_normal_iterations_min = 1;
const num_t bilateralnormal_num_normal_iterations_max = 100;
const num_t bilateralnormal_num_normal_iterations_delta = 1;

const num_t bilateralnormal_num_vertex_iterations_default = 10;
const num_t bilateralnormal_num_vertex_iterations_min = 1;
const num_t bilateralnormal_num_vertex_iterations_max = 100;
const num_t bilateralnormal_num_vertex_iterations_delta = 1;

const num_t bilateralnormal_sigma_c_scalar_default = 1.0;
const num_t bilateralnormal_sigma_c_scalar_min = 0.5;
const num_t bilateralnormal_sigma_c_scalar_max = 5.0;
const num_t bilateralnormal_sigma_c_scalar_delta = 0.01;

const num_t bilateralnormal_sigma_s_default = 0.4;
const num_t bilateralnormal_sigma_s_min = 0.1;
const num_t bilateralnormal_sigma_s_max = 0.9;
const num_t bilateralnormal_sigma_s_delta = 0.01;

const string bilateralnormal_parameter_order = "<number of normal iterations> <number of vertex iterations> <sigma_c scalar> <sigma_s>";

const vector<num_t>
    bilateralnormal_parameters_default = {bilateralnormal_num_normal_iterations_default,bilateralnormal_num_vertex_iterations_default,bilateralnormal_sigma_c_scalar_default,bilateralnormal_sigma_s_default};
const vector<num_t>
    bilateralnormal_parameters_min = {bilateralnormal_num_normal_iterations_min,bilateralnormal_num_vertex_iterations_min,bilateralnormal_sigma_c_scalar_min,bilateralnormal_sigma_s_min};
const vector<num_t>
    bilateralnormal_parameters_max = {bilateralnormal_num_normal_iterations_max,bilateralnormal_num_vertex_iterations_max,bilateralnormal_sigma_c_scalar_max,bilateralnormal_sigma_s_max};
const vector<num_t>
    bilateralnormal_parameters_delta = {bilateralnormal_num_normal_iterations_delta,bilateralnormal_num_vertex_iterations_delta,bilateralnormal_sigma_c_scalar_delta,bilateralnormal_sigma_s_delta};

//Guided
const string guided_key_word = "guided";

const num_t guided_num_normal_iterations_default = 20;
const num_t guided_num_normal_iterations_min = 1;
const num_t guided_num_normal_iterations_max = 100;
const num_t guided_num_normal_iterations_delta = 1;

const num_t guided_num_vertex_iterations_default = 10;
const num_t guided_num_vertex_iterations_min = 1;
const num_t guided_num_vertex_iterations_max = 100;
const num_t guided_num_vertex_iterations_delta = 1;

const num_t guided_sigma_c_scalar_default = 1.0;
const num_t guided_sigma_c_scalar_min = 0.5;
const num_t guided_sigma_c_scalar_max = 5.0;
const num_t guided_sigma_c_scalar_delta = 0.01;

const num_t guided_sigma_s_default = 0.4;
const num_t guided_sigma_s_min = 0.1;
const num_t guided_sigma_s_max = 0.9;
const num_t guided_sigma_s_delta = 0.01;

const num_t guided_radius_scalar_default = 2.0;
const num_t guided_radius_scalar_min = 1.0;
const num_t guided_radius_scalar_max = 10.0;
const num_t guided_radius_scalar_delta = 0.01;

const string guided_parameter_order = "<number of normal iterations> <number of vertex iterations> <sigma_c scalar> <sigma_s> <radius scalar>";

const vector<num_t> guided_parameters_default =
{guided_num_normal_iterations_default,guided_num_vertex_iterations_default,guided_sigma_c_scalar_default,guided_sigma_s_default,guided_radius_scalar_default};
const vector<num_t> guided_parameters_min =
{guided_num_normal_iterations_min,guided_num_vertex_iterations_min,guided_sigma_c_scalar_min,guided_sigma_s_min,guided_radius_scalar_min};
const vector<num_t> guided_parameters_max =
{guided_num_normal_iterations_max,guided_num_vertex_iterations_max,guided_sigma_c_scalar_max,guided_sigma_s_max,guided_radius_scalar_max};
const vector<num_t> guided_parameters_delta =
{guided_num_normal_iterations_delta,guided_num_vertex_iterations_delta,guided_sigma_c_scalar_delta,guided_sigma_s_delta,guided_radius_scalar_delta};

inline void printKeyWords()
{
    printf("\nalgorithms: %s | %s | %s | %s | %s | %s | %s \n", uniform_key_word.c_str(), taubin_key_word.c_str(), hc_key_word.c_str(), bilateral_key_word.c_str(), fast_key_word.c_str(), bilateralnormal_key_word.c_str(), guided_key_word.c_str());
}

inline void printUsage()
{
    printf("\nusage:\n<program_name> <input_file_name> <output_file_name> <denoising_algorithm> <parameter_1> ... <parameter_n> \n");
    printf(" Algorithm : parameters \n");
    printf("    %s : %s \n",uniform_key_word.c_str(),uniform_parameter_order.c_str());
    printf("    %s : %s \n",taubin_key_word.c_str(),taubin_parameter_order.c_str());
    printf("    %s : %s \n",hc_key_word.c_str(),hc_parameter_order.c_str());
    printf("    %s : %s \n",bilateral_key_word.c_str(),bilateral_parameter_order.c_str());
    printf("    %s : %s \n",fast_key_word.c_str(),fast_parameter_order.c_str());
    printf("    %s : %s \n",bilateralnormal_key_word.c_str(),bilateralnormal_parameter_order.c_str());
    printf("    %s : %s \n",guided_key_word.c_str(),guided_parameter_order.c_str());
}

#endif // CUSTOM_H
