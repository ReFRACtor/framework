#include "error_analysis.h"
#include "absorber_absco.h"
#include "unit_test_support.h"
#include "solver_finished_fixture.h"

using namespace FullPhysics;
using namespace blitz;
BOOST_FIXTURE_TEST_SUITE(error_analysis, SolverFinishedFixture)

BOOST_AUTO_TEST_CASE(basic)
{
    if (false) {
        // Enable to output expected values
        std::cerr << std::setprecision(20) << std::scientific;
        std::cerr << "double signal_level_0 = " << error_analysis->signal_level(0)<< ";" << std::endl;
        std::cerr << "double signal_level_1 = " << error_analysis->signal_level(1)<< ";" << std::endl;
        std::cerr << "double signal_level_2 = " << error_analysis->signal_level(2)<< ";" << std::endl;
        std::cerr << "double noise_level_0 = " << error_analysis->noise_level(0)<< ";" << std::endl;
        std::cerr << "double noise_level_1 = " << error_analysis->noise_level(1)<< ";" << std::endl;
        std::cerr << "double noise_level_2 = " << error_analysis->noise_level(2)<< ";" << std::endl;
        std::cerr << "double residual_sum_sq_0 = " << error_analysis->residual_sum_sq(0)<< ";" << std::endl;
        std::cerr << "double residual_sum_sq_1 = " << error_analysis->residual_sum_sq(1)<< ";" << std::endl;
        std::cerr << "double residual_sum_sq_2 = " << error_analysis->residual_sum_sq(2)<< ";" << std::endl;
        std::cerr << "double residual_mean_sq_0 = " << error_analysis->residual_mean_sq(0)<< ";" << std::endl;
        std::cerr << "double residual_mean_sq_1 = " << error_analysis->residual_mean_sq(1)<< ";" << std::endl;
        std::cerr << "double residual_mean_sq_2 = " << error_analysis->residual_mean_sq(2)<< ";" << std::endl;
        std::cerr << "double reduced_chisq_0 = " << error_analysis->reduced_chisq(0)<< ";" << std::endl;
        std::cerr << "double reduced_chisq_2 = " << error_analysis->reduced_chisq(2)<< ";" << std::endl;
        std::cerr << "double reduced_chisq_1 = " << error_analysis->reduced_chisq(1)<< ";" << std::endl;
        std::cerr << "double degrees_of_freedom = " << error_analysis->degrees_of_freedom()<< ";" << std::endl;
    }

    double signal_level_0 = 6.30019854852392222720e+19;
    double signal_level_1 = 2.42517250040087429120e+19;
    double signal_level_2 = 1.06198430635953991680e+19;
    double noise_level_0 = 1.95667167800942048000e+17;
    double noise_level_1 = 5.24157317853121760000e+16;
    double noise_level_2 = 3.75894412497323760000e+16;
    double residual_sum_sq_0 = 3.67542735684385317887e+37;
    double residual_sum_sq_1 = 1.47588157224684958536e+37;
    double residual_sum_sq_2 = 3.02456794443312232363e+36;
    double residual_mean_sq_0 = 2.09676994607205792000e+17;
    double residual_mean_sq_1 = 1.50338137113825120000e+17;
    double residual_mean_sq_2 = 6.10314457471845440000e+16;
    double reduced_chisq_0 = 4.09113394101837091199e+00;
    double reduced_chisq_2 = 4.06877095551921019734e+00;
    double reduced_chisq_1 = 7.98982656306322880368e+00;
    double degrees_of_freedom = 1.84488824435404410451e+01;

    BOOST_CHECK_CLOSE(error_analysis->signal_level(0), signal_level_0, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->signal_level(1), signal_level_1, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->signal_level(2), signal_level_2, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->noise_level(0), noise_level_0, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->noise_level(1), noise_level_1, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->noise_level(2), noise_level_2, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->residual_sum_sq(0), residual_sum_sq_0, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->residual_sum_sq(1), residual_sum_sq_1, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->residual_sum_sq(2), residual_sum_sq_2, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->residual_mean_sq(0), residual_mean_sq_0, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->residual_mean_sq(1), residual_mean_sq_1, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->residual_mean_sq(2), residual_mean_sq_2, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->relative_residual_mean_sq(0), error_analysis->residual_mean_sq(0) / error_analysis->signal_level(0), 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->relative_residual_mean_sq(1), error_analysis->residual_mean_sq(1) / error_analysis->signal_level(1), 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->relative_residual_mean_sq(2), error_analysis->residual_mean_sq(2) / error_analysis->signal_level(2), 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->reduced_chisq(0), reduced_chisq_0, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->reduced_chisq(2), reduced_chisq_2, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->reduced_chisq(1), reduced_chisq_1, 1e-8);
    BOOST_CHECK_CLOSE(error_analysis->degrees_of_freedom(), degrees_of_freedom, 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
