#include "solver_iteration_log.h"

using namespace FullPhysics;
using namespace blitz;

void SolverIterationLog::notify_update(const IterativeSolver& solver)
{ 
    // Iteration index 0 is the initial state
    int iter_index = solver.num_accepted_steps();

    if (iter_index < 0) {
        throw Exception("Number of iteration indexes indicated by IterativeSolver is negative");
    }

    blitz::Array<std::string, 1> sv_names = sv_obj->state_vector_name();

    blitz::Array<double, 1> curr_sv(solver.accepted_points().at(iter_index));
    blitz::Array<double, 1> dx;

    if (iter_index > 0) {
        dx.resize(curr_sv.rows());
        dx = solver.accepted_points().at(iter_index) - solver.accepted_points().at(iter_index - 1);
    }

    std::stringstream header_log;
    header_log << std::endl;
    if (iter_index > 0) {
        header_log << "Accepted iteration #" << iter_index << std::endl;
    } else {
        header_log << "Initial retrieval state" << std::endl;
    }
    header_log << "Solver status: " << solver.status_str() << std::endl;
    header_log << "Cost function value: " << solver.cost_at_accepted_points().at(iter_index) << std::endl;
    Logger::info() << header_log.str();

    std::stringstream sv_log;
    sv_log << std::endl
           << std::setw(SV_PRINT_WIDTH) 
           << "State Vector";
    if (iter_index > 0) {
        sv_log << std::setw(SV_PRINT_WIDTH) 
               << "Dx";
    }
    sv_log << "  Name" << std::endl;
    for(int sv_idx = 0; sv_idx < curr_sv.rows(); sv_idx++) {
        sv_log << std::setprecision(SV_PRINT_WIDTH-7) // leave room for -/+, ., exponent and space, etc
               << std::setw(SV_PRINT_WIDTH) 
               << curr_sv(sv_idx);
        if (iter_index > 0) {
            sv_log << std::setw(SV_PRINT_WIDTH) 
                   << dx(sv_idx);
        }
        sv_log << "  "
               << sv_names(sv_idx)
               << std::endl;
    }
    Logger::info() << sv_log.str();
}
