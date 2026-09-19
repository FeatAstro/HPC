#include <H5Cpp.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using ScalarFunction = std::function<double(double)>;

double derivative_first_order(const ScalarFunction& f, double x, double h) {
    return (f(x + h) - f(x)) / h;
}

double derivative_second_order(const ScalarFunction& f, double x, double h) {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}

// n_points values spaced evenly in log10 between 10^log_start and 10^log_end,
// so that error-vs-h behaves as a straight line on a log-log plot.
std::vector<double> log_space(double log_start, double log_end, int n_points) {
    std::vector<double> values(n_points);
    for (int i = 0; i < n_points; ++i) {
        double t = static_cast<double>(i) / (n_points - 1);
        double log_value = log_start + t * (log_end - log_start);
        values[i] = std::pow(10.0, log_value);
    }
    return values;
}

struct ConvergenceResults {
    std::vector<double> h;
    std::vector<double> error_first_order;
    std::vector<double> error_second_order;
};

// Evaluates the absolute error of both finite-difference schemes against the
// exact derivative, for every step size in h_values.
ConvergenceResults compute_convergence(const ScalarFunction& f, double exact_derivative,
                                        double x, const std::vector<double>& h_values) {
    ConvergenceResults results;
    results.h = h_values;
    results.error_first_order.resize(h_values.size());
    results.error_second_order.resize(h_values.size());

    for (std::size_t i = 0; i < h_values.size(); ++i) {
        double h = h_values[i];
        double approx_first_order = derivative_first_order(f, x, h);
        double approx_second_order = derivative_second_order(f, x, h);

        results.error_first_order[i] = std::abs(approx_first_order - exact_derivative);
        results.error_second_order[i] = std::abs(approx_second_order - exact_derivative);
    }

    return results;
}

void write_convergence_results(const std::string& filename, const ConvergenceResults& results,
                                double x, double exact_derivative) {
    H5::H5File file(filename, H5F_ACC_TRUNC);

    hsize_t dims[1] = {static_cast<hsize_t>(results.h.size())};
    H5::DataSpace dataspace(1, dims);

    auto write_dataset = [&](const std::string& name, const std::vector<double>& data) {
        H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(data.data(), H5::PredType::NATIVE_DOUBLE);
    };

    write_dataset("h", results.h);
    write_dataset("error_first_order", results.error_first_order);
    write_dataset("error_second_order", results.error_second_order);

    H5::DataSpace scalar_space(H5S_SCALAR);
    H5::Attribute attr_x = file.createAttribute("x", H5::PredType::NATIVE_DOUBLE, scalar_space);
    attr_x.write(H5::PredType::NATIVE_DOUBLE, &x);
    H5::Attribute attr_exact = file.createAttribute("exact_derivative", H5::PredType::NATIVE_DOUBLE, scalar_space);
    attr_exact.write(H5::PredType::NATIVE_DOUBLE, &exact_derivative);
}

int main() {
    // sin/cos is used because its derivative is known analytically, so the
    // finite-difference approximations can be checked against an exact value.
    ScalarFunction f = [](double x) { return std::sin(x); };
    const double x = 1.0;
    const double exact_derivative = std::cos(x);

    // h spans 10^-1 down to 10^-12: wide enough to see both the truncation-error
    // slope (large h) and the round-off floor where it plateaus (small h).
    const int n_points = 56;
    const double log_h_start = -1.0;
    const double log_h_end = -12.0;
    std::vector<double> h_values = log_space(log_h_start, log_h_end, n_points);

    ConvergenceResults results = compute_convergence(f, exact_derivative, x, h_values);

    std::cout << "Computed " << n_points << " finite-difference errors for h in [1e"
              << log_h_end << ", 1e" << log_h_start << "]" << std::endl;

    const std::string filename = "results.h5";
    write_convergence_results(filename, results, x, exact_derivative);

    std::cout << "Results written to " << filename << std::endl;

    return 0;
}
