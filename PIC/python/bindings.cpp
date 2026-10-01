// Exposes the C++ core to Python as the module `pic`.
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "boris.hpp"
#include "field_solver.hpp"
#include "grid.hpp"
#include "loading.hpp"
#include "population.hpp"
#include "vec3.hpp"

namespace py = pybind11;

static py::array_t<double> to_numpy(const std::vector<double>& values) {
    return py::array_t<double>(values.size(), values.data());
}

// forcecast: any numeric numpy array is converted to contiguous doubles
using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
static void from_numpy(std::vector<double>& destination, const DoubleArray& array) {
    destination.assign(array.data(), array.data() + array.size());
}

static std::vector<double> to_vector(const DoubleArray& array) {
    return std::vector<double>(array.data(), array.data() + array.size());
}

static std::string vec3_repr(const Vec3& a) {
    return "Vec3(" + std::to_string(a.x) + ", " + std::to_string(a.y) + ", " + std::to_string(a.z) + ")";
}

PYBIND11_MODULE(pic, m) {
    m.doc() = "1D3V particle-in-cell core (C++), wrapped for testing";

    py::class_<Vec3>(m, "Vec3")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z)
        .def("__repr__", &vec3_repr);

    // Step 1: one particle, Boris pusher
    py::class_<Particle>(m, "Particle")
        .def(py::init<>())
        .def_readwrite("x", &Particle::x)
        .def_readwrite("v", &Particle::v)
        .def_readwrite("q", &Particle::q)
        .def_readwrite("m", &Particle::m);

    py::class_<Trajectory>(m, "Trajectory")
        .def_property_readonly("t", [](const Trajectory& traj) { return to_numpy(traj.t); })
        .def_property_readonly("x", [](const Trajectory& traj) { return to_numpy(traj.x); })
        .def_property_readonly("vx", [](const Trajectory& traj) { return to_numpy(traj.vx); })
        .def_property_readonly("vy", [](const Trajectory& traj) { return to_numpy(traj.vy); })
        .def_property_readonly("vz", [](const Trajectory& traj) { return to_numpy(traj.vz); });

    m.def("run_particle", &run_particle,
          py::arg("particle"), py::arg("E"), py::arg("B"), py::arg("dt"), py::arg("n_steps"));

    // Steps 2b and 3a: periodic grid, Yee layout, order-1 deposit and gather.
    py::enum_<Centering>(m, "Centering")
        .value("node", Centering::node)
        .value("cell_centre", Centering::cell_centre);

    py::class_<Grid>(m, "Grid")
        .def(py::init<int, double>(), py::arg("nx"), py::arg("dx"))
        .def_readonly("nx", &Grid::nx)
        .def_readonly("dx", &Grid::dx)
        .def("length", &Grid::length)
        .def("node_position", &Grid::node_position, py::arg("i"))
        .def("cell_centre_position", &Grid::cell_centre_position, py::arg("i"));

    py::class_<VectorField>(m, "VectorField")
        .def(py::init<const Grid&>(), py::arg("grid"))
        .def(py::init<const Grid&, Centering, Centering, Centering>(), py::arg("grid"),
             py::arg("x_centering"), py::arg("y_centering"), py::arg("z_centering"))
        .def_readonly("x_centering", &VectorField::x_centering)
        .def_readonly("y_centering", &VectorField::y_centering)
        .def_readonly("z_centering", &VectorField::z_centering)
        .def_property("x", [](const VectorField& field) { return to_numpy(field.x); },
                      [](VectorField& field, const DoubleArray& values) { from_numpy(field.x, values); })
        .def_property("y", [](const VectorField& field) { return to_numpy(field.y); },
                      [](VectorField& field, const DoubleArray& values) { from_numpy(field.y, values); })
        .def_property("z", [](const VectorField& field) { return to_numpy(field.z); },
                      [](VectorField& field, const DoubleArray& values) { from_numpy(field.z, values); });

    py::class_<Fields>(m, "Fields")
        .def(py::init<const Grid&>(), py::arg("grid"))
        .def_readwrite("E", &Fields::E)
        .def_readwrite("B", &Fields::B);

    m.def("gather", [](const VectorField& field, const Grid& grid, double x) {
        return gather(field, particle_shape(grid, x));
    }, py::arg("field"), py::arg("grid"), py::arg("x"));

    m.def("derivative_at_cell_centres", [](const DoubleArray& node_values, const Grid& grid) {
        return to_numpy(derivative_at_cell_centres(to_vector(node_values), grid));
    }, py::arg("node_values"), py::arg("grid"));
    m.def("derivative_at_nodes", [](const DoubleArray& cell_centre_values, const Grid& grid) {
        return to_numpy(derivative_at_nodes(to_vector(cell_centre_values), grid));
    }, py::arg("cell_centre_values"), py::arg("grid"));
    m.def("average_at_cell_centres", [](const DoubleArray& node_values) {
        return to_numpy(average_at_cell_centres(to_vector(node_values)));
    }, py::arg("node_values"));
    m.def("average_at_nodes", [](const DoubleArray& cell_centre_values) {
        return to_numpy(average_at_nodes(to_vector(cell_centre_values)));
    }, py::arg("cell_centre_values"));

    // Steps 2a-2b: N particles. Reading pop.x returns a copy, `pop.x = array` copies into C++.
    py::class_<Population>(m, "Population")
        .def(py::init<const std::string&, double, double>(), py::arg("name"), py::arg("q"), py::arg("m"))
        .def_readonly("name", &Population::name)
        .def_readwrite("q", &Population::q)
        .def_readwrite("m", &Population::m)
        .def_property("x", [](const Population& pop) { return to_numpy(pop.x); },
                      [](Population& pop, const DoubleArray& values) { from_numpy(pop.x, values); })
        .def_property("vx", [](const Population& pop) { return to_numpy(pop.vx); },
                      [](Population& pop, const DoubleArray& values) { from_numpy(pop.vx, values); })
        .def_property("vy", [](const Population& pop) { return to_numpy(pop.vy); },
                      [](Population& pop, const DoubleArray& values) { from_numpy(pop.vy, values); })
        .def_property("vz", [](const Population& pop) { return to_numpy(pop.vz); },
                      [](Population& pop, const DoubleArray& values) { from_numpy(pop.vz, values); })
        .def_property("w", [](const Population& pop) { return to_numpy(pop.w); },
                      [](Population& pop, const DoubleArray& values) { from_numpy(pop.w, values); })
        .def_property_readonly("density", [](const Population& pop) { return to_numpy(pop.density); })
        .def_readonly("flux", &Population::flux)
        .def_property_readonly("kinetic_energy_density",
                               [](const Population& pop) { return to_numpy(pop.kinetic_energy_density); })
        .def("__len__", &Population::size)
        .def("push", &Population::push, py::arg("grid"), py::arg("fields"), py::arg("dt"))
        .def("deposit", &Population::deposit, py::arg("grid"))
        .def("temperature", [](const Population& pop) { return to_numpy(pop.temperature()); });

    m.def("total_density", [](const std::vector<const Population*>& populations, const Grid& grid) {
        return to_numpy(total_density(populations, grid));
    }, py::arg("populations"), py::arg("grid"));
    m.def("bulk_velocity", &bulk_velocity, py::arg("populations"), py::arg("grid"));

    // Step 2c: loading from user profiles (Python functions of x).
    m.def("load_maxwellian", &load_maxwellian, py::arg("population"), py::arg("grid"), py::arg("particles_per_cell"),
          py::arg("density"), py::arg("bulk_velocity"), py::arg("temperature"), py::arg("seed"));

    py::class_<PopulationHistory>(m, "PopulationHistory")
        .def_property_readonly("t", [](const PopulationHistory& history) { return to_numpy(history.t); })
        .def_property_readonly("kinetic_energy",
                               [](const PopulationHistory& history) { return to_numpy(history.kinetic_energy); })
        .def_property_readonly("mean_vx", [](const PopulationHistory& history) { return to_numpy(history.mean_vx); })
        .def_property_readonly("mean_vy", [](const PopulationHistory& history) { return to_numpy(history.mean_vy); })
        .def_property_readonly("mean_vz", [](const PopulationHistory& history) { return to_numpy(history.mean_vz); });

    m.def("run_population", &run_population, py::arg("population"), py::arg("grid"), py::arg("fields"),
          py::arg("dt"), py::arg("n_steps"));

    // Step 3b: Ampere and Ohm.
    py::class_<PlasmaParameters>(m, "PlasmaParameters")
        .def(py::init<>())
        .def_readwrite("electron_temperature", &PlasmaParameters::electron_temperature)
        .def_readwrite("resistivity", &PlasmaParameters::resistivity)
        .def_readwrite("hyper_resistivity", &PlasmaParameters::hyper_resistivity);

    m.def("ampere", &ampere, py::arg("B"), py::arg("grid"));
    m.def("ohm", [](const DoubleArray& density, const VectorField& bulk_velocity, const VectorField& B,
                    const VectorField& j, const Grid& grid, const PlasmaParameters& parameters) {
        return ohm(to_vector(density), bulk_velocity, B, j, grid, parameters);
    }, py::arg("density"), py::arg("bulk_velocity"), py::arg("B"), py::arg("j"), py::arg("grid"),
       py::arg("parameters"));

    // Step 3c: Faraday and the iterated Crank-Nicolson time step, ions held fixed.
    m.def("advance_fields_with_fixed_ions",
          [](Fields& fields, const DoubleArray& density, const VectorField& bulk_velocity, const Grid& grid,
             const PlasmaParameters& parameters, double dt, int n_steps) {
              const std::vector<double> ion_density = to_vector(density);
              for (int step = 0; step < n_steps; ++step) {
                  advance_fields_with_fixed_ions(fields, ion_density, bulk_velocity, grid, parameters, dt);
              }
          },
          py::arg("fields"), py::arg("density"), py::arg("bulk_velocity"), py::arg("grid"), py::arg("parameters"),
          py::arg("dt"), py::arg("n_steps") = 1);
}
