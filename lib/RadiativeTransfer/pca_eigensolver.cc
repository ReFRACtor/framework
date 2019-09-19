// PCA using Eigenproblem methods

#include "pca_eigensolver.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void prepare_eigenmatrix(int *m, int *n, int *t, double *x, double* y, double *ccm);
    void pca_asymtx(double *aad, int *m, int *ia, int *ievec, int *ievec2, double *tol, double *evecd, double *evald, int *ier, double *wkd, int *message_len, char *message, bool *bad_status);
    void pca_ranker(int *n, double* arrin, int *indx);

    void pca_eigensolver(int *Max_Eofs, int *maxpoints, int *maxlayers, int *maxlayers2, int *n_Eofs, int *npoints, int *nlayers, int *nlayers2, double *taudp, double *omega, double *Atmosmean, double *Eofs, double *PrinComps, bool *fail, int *message_len, char *message, int *trace_len, char *trace);
    void pca_eigensolver_alb(int *Max_Eofs, int *maxpoints, int *maxlayers, int *maxlayers21, int *n_Eofs, int *npoints, int *nlayers, int *nlayers2, int *nlayers21, double *taudp, double *omega, double *albedo, double *Atmosmean, double *Albmean, double *Eofs, double *PrinComps, bool *fail, int *message_len, char *message, int *trace_len, char *trace);

}

std::vector<blitz::Array<double, 3> > PCAEigenSolver::data_perturbations() const
{
    auto all = Range::all();
    int num_vars = data_mean().size();
    int num_eofs = eof_properties()[0].cols();

    std::vector<Array<double, 3> > result;
    for(int ivar = 0; ivar < num_vars; ivar++) {
        int num_data = data_mean()[ivar].rows();
        Array<double, 3> var_pert(num_data, 2, num_eofs);
        Array<double, 1> var_mean(data_mean()[ivar]);
        Array<double, 2> var_eof(eof_properties()[ivar]);

        for(int eof_idx = 0; eof_idx < num_eofs; eof_idx++) {
            var_pert(all, 0, eof_idx) = var_mean + var_eof(all, eof_idx);
            var_pert(all, 1, eof_idx) = var_mean - var_eof(all, eof_idx);
        }

        result.push_back(var_pert);
    }

    return result;
}

blitz::Array<double, 2> PCAEigenSolver::correction(const blitz::Array<double, 1>& lidort_mean, 
                                                   const blitz::Array<double, 1>& twostream_mean, 
                                                   const blitz::Array<double, 1>& first_order_mean,
                                                   const blitz::Array<double, 2>& lidort_plus,
                                                   const blitz::Array<double, 2>& twostream_plus,
                                                   const blitz::Array<double, 2>& first_order_plus,
                                                   const blitz::Array<double, 2>& lidort_minus,
                                                   const blitz::Array<double, 2>& twostream_minus,
                                                   const blitz::Array<double, 2>& first_order_minus)
{
    Range all = Range::all();
    int num_eofs = principal_components().rows();
    int num_points = principal_components().cols();
    int num_stokes = lidort_mean.rows();

    // Check that plus/minus have appropriate shape
    if(lidort_plus.rows() != num_eofs || lidort_minus.rows() != num_eofs) {
        throw Exception("LIDORT radiances values for EOF correction do not match number of EOFs");
    }

    if(twostream_plus.rows() != num_eofs || twostream_minus.rows() != num_eofs) {
        throw Exception("2stream radiances values for EOF correction do not match number of EOFs");
    }

    if(first_order_plus.rows() != num_eofs || first_order_minus.rows() != num_eofs) {
        throw Exception("First Order radiances values for EOF correction do not match number of EOFs");
    }

    if(lidort_mean.rows() != num_stokes || lidort_plus.cols() != num_stokes || lidort_minus.cols() != num_stokes) {
        throw Exception("LIDORT number of stokes is inconsistent");
    }

    if(twostream_mean.rows() != num_stokes || twostream_plus.cols() != num_stokes || twostream_minus.cols() != num_stokes) {
        throw Exception("2stream number of stokes is inconsistent");
    }

    if(first_order_mean.rows() != num_stokes || first_order_plus.cols() != num_stokes || first_order_minus.cols() != num_stokes) {
        throw Exception("First Order number of stokes is inconsistent");
    }

    firstIndex i1;
    Array<double, 1> hi_ss( lidort_mean + first_order_mean );
    Array<double, 1> lo_ss( twostream_mean + first_order_mean );
    Array<double, 1> id_mean( where(lo_ss(i1) > 0, log( hi_ss / lo_ss ), 0.0) );

    Array<double, 2> term_1(num_eofs, num_stokes);
    Array<double, 2> term_2(num_eofs, num_stokes);

    for(int eof_idx = 0; eof_idx < num_eofs; eof_idx++) {
        Array<double, 1> hi_ss_plus( lidort_plus(eof_idx, all) + first_order_plus(eof_idx, all) );
        Array<double, 1> lo_ss_plus( twostream_plus(eof_idx, all) + first_order_plus(eof_idx, all) );
        Array<double, 1> id_plus( where(lo_ss_plus(i1) > 0, log(hi_ss_plus / lo_ss_plus), 0.0) );

        Array<double, 1> hi_ss_minus( lidort_minus(eof_idx, all) + first_order_minus(eof_idx, all) );
        Array<double, 1> lo_ss_minus( twostream_minus(eof_idx, all) + first_order_minus(eof_idx, all) );
        Array<double, 1> id_minus( where(lo_ss_minus(i1) > 0, log(hi_ss_minus / lo_ss_minus), 0.0) );
  
        term_1(eof_idx, all) = (id_plus - id_minus) / 2;
        term_2(eof_idx, all) = (id_plus + id_minus - 2*id_mean) / 2;
    }

    Array<double, 2> correction(num_points, num_stokes);
    for(int point_idx = 0; point_idx < num_points; point_idx++) {
        Array<double, 1> ieof(num_stokes);
        ieof = id_mean;

        for(int eof_idx = 0; eof_idx < num_eofs; eof_idx++) {
            double x = principal_components()(eof_idx, point_idx);
            ieof = ieof + term_1(eof_idx, all) * x + term_2(eof_idx, all) * x * x;
        }
        correction(point_idx, all) = exp(ieof);
    }

    return correction;
}

// ----------

PCAEigenSolverGeneric::PCAEigenSolverGeneric(const std::vector<Array<double, 2> >& gridded_data, int num_eofs)
{
    solve(gridded_data, num_eofs);
}

// Adapted from the original Fortran code to be more generic

void PCAEigenSolverGeneric::solve(const std::vector<Array<double, 2> >& gridded_data, int num_eofs)
{
    Range all = Range::all();
    secondIndex i2;

    int num_vars = gridded_data.size();

    // Figure out number of data values to pack from the input data
    // Expected that the second dimension of each of the variable should
    // be the number of grid points and all variables should be the same
    int num_packed = 0;
    int num_points = 0;
    for(int ivar = 0; ivar < num_vars; ivar++) {
        // Ensure that all variables use the same number of points
        if(num_points == 0) {
            num_points = gridded_data[ivar].cols();
        } else if (num_points != gridded_data[ivar].cols()) {
            Exception err;
            err << "The variable at index " << ivar << " "
                << "does not have the same number of points " << num_points << " "
                << "as other variables, instead it has " << gridded_data[ivar].cols() << " points";
            throw err;
        }

        num_packed += gridded_data[ivar].rows();
    }

    // 1. Get & process data
    // ==================
    Array<double, 2> packed_data(num_packed, num_points, ColumnMajorArray<2>());

    // process data into one array, take logarithm
    int ipacked = 0;
    for(int ivar = 0; ivar < num_vars; ivar++) {
        int num_data = gridded_data[ivar].rows();
        Range packed_r(ipacked, ipacked + num_data - 1);
        packed_data(packed_r, all) = log(gridded_data[ivar](all, all));
        ipacked += num_data;
    }

    // Remove time-mean for each packed variable, mean over data index
    for(int ipack = 0; ipack < num_packed; ipack++) {
        packed_data(ipack, all) = packed_data(ipack, all) - mean(packed_data(ipack, all));
    }

    // Compute log mean value over the spectral grid for each varible type
    for(int ivar = 0; ivar < num_vars; ivar++) {
        int num_data = gridded_data[ivar].rows();

        Array<double, 1> var_mean(num_data);
        var_mean = mean(log(gridded_data[ivar]), i2);
        atmos_mean_.push_back(var_mean);
    }

    // 2. Prepare Eigenmatrix and Solve Eigenproblem
    // ==========================================

    // Prepare matrix (cross covariances)
    Array<double, 2> eigenmat(num_packed, num_packed, ColumnMajorArray<2>());
    Array<double, 2> eigenmat_save(num_packed, num_packed, ColumnMajorArray<2>());

    prepare_eigenmatrix(&num_packed, &num_packed, &num_points, packed_data.dataFirst(), packed_data.dataFirst(), eigenmat.dataFirst());

    // Save eigenmat in case tol needs to be changed
    eigenmat_save = eigenmat;

    // tolerance input to Eigenpackage module ASYMTX
    double tol = 1.0e-6;

    // (output from Eigenpackage module ASYMTX)
    Array<double, 2> evec(num_packed, num_packed, ColumnMajorArray<2>());
    Array<double, 1> ksq(num_packed);
    Array<double, 1> wk(2 * num_packed);
    bool asymtx_failure = false;
    
    int message_len = 100;
    blitz::Array<char, 1> asymtx_message(message_len);

    // Solve eigenproblem using PCA_ASYMTX
    int ier = 999;
    int num_packed2 = 2*num_packed;
    while(ier > 0 && not asymtx_failure) {
        pca_asymtx(eigenmat.dataFirst(), &num_packed, &num_packed, &num_packed, &num_packed2, &tol,
               evec.dataFirst(), ksq.dataFirst(), &ier, wk.dataFirst(), &message_len, asymtx_message.dataFirst(), &asymtx_failure);

        // Change tolerance and rerun PCA_ASYMTX if eigenvalue has not converged
        if(ier > 0) {
            tol = tol * 10.0;
            eigenmat = eigenmat_save;
        }
    }

    // Exception handling 1
    if (asymtx_failure) {
        Exception err;
        std::string message_str;
        message_str = std::string(asymtx_message(all).begin(), asymtx_message(all).end());
        err << "pca_asymtx failed with message: " << message_str;
        throw err;
    }

    // Exception handling 2
    if (ier > 0) {
        Exception err;
        err << "pca_asymtx error: eigenvalue " << ier << " has not converged";
        throw err;
    }

    // Normalize vectors
    for(int ipack2 = 0; ipack2 < num_packed; ipack2++) {
        double norm = sum( evec(all, ipack2) * evec(all, ipack2) );
        for(int ipack1 = 0; ipack1 < num_packed; ipack1++) {
            evec(ipack1, ipack2) = evec(ipack1, ipack2) / norm;
        }
    }

    // rank absolute(eigenvalues). [subroutine "Ranker" is same as Indexx1]
    Array<double, 1> ksq_abs(num_packed);
    ksq_abs = abs(ksq);

    Array<int, 1> order(num_packed);

    pca_ranker(&num_packed, ksq_abs.dataFirst(), order.dataFirst());

    // Rank the Eigenvectors
    Array<double, 1> ksq_ordered(num_packed);
    Array<double, 2> evec_2(num_packed, num_packed, ColumnMajorArray<2>());

    for(int ipack2 = 0; ipack2 < num_packed; ipack2++) {
        int iord = num_packed - ipack2 - 1;
        ksq_ordered(ipack2) = ksq_abs(order(iord)-1);

        for(int ipack1 = 0; ipack1 < num_packed; ipack1++) {
            evec_2(ipack1, ipack2)= evec(ipack1, order(iord)-1);
        }
    }

    // 3. Prepare EOF and PC output 
    // =========================
    
    // Resize EOFS and Principal Components
    Array<double, 2> eofs_packed(num_eofs, num_packed);
    prin_comps_.resize(num_eofs, num_points);

    // *** Only perform for the first few EOFs
    for(int aa = 0; aa < num_eofs; aa++) {

        // set the usable eigenvalue
        double lambda = ksq_ordered(aa);

        // Normalize everything according to SQRT(Lambda)
        double stdv = std::sqrt(lambda);

        // EOFs (Unnormalized) --> Transpose the Eigenvectors
        eofs_packed(aa, all) = evec_2(all, aa);

        // Project data onto E1 basis
        // -- Set the principal components (unnormalized)
        for(int w = 0; w < num_points; w++) {
            prin_comps_(aa, w) = sum(eofs_packed(aa, all) * packed_data(all,w));
        }

        // Final normalization of EOFs and PCs
        eofs_packed(aa, all) = eofs_packed(aa, all) * stdv;
        prin_comps_(aa, all) = prin_comps_(aa, all) / stdv;
    }

    // Unpack eofs array values
    ipacked = 0;
    for(int ivar = 0; ivar < num_vars; ivar++) {
        int num_data = gridded_data[ivar].rows();
        Range packed_r(ipacked, ipacked + num_data - 1);
    
        Array<double, 2> data_eofs(num_data, num_eofs);
        data_eofs = eofs_packed(all, packed_r).transpose(secondDim, firstDim);
        eofs_.push_back( data_eofs );
        ipacked += num_data;
    }

}

// ----------

PCAEigenSolverFortran::PCAEigenSolverFortran(const blitz::Array<double, 3>& gridded_data, int num_eofs)
{
    solve(gridded_data, num_eofs);
}
    
std::vector<blitz::Array<double, 1> > PCAEigenSolverFortran::data_mean() const
{
    std::vector<blitz::Array<double, 1> > result;
    for (int ivar = 0; ivar < atmos_mean_.cols(); ivar++) {
        result.push_back( Array<double, 1>(atmos_mean_(Range::all(), ivar)) );
    }
    return result;
}

std::vector<blitz::Array<double, 2> > PCAEigenSolverFortran::eof_properties() const
{
    auto all = Range::all();
    int num_layers = atmos_mean_.rows();
    int num_eofs = eofs_.rows();

    int ipacked = 0;
    std::vector<blitz::Array<double, 2> > result;
    for (int ivar = 0; ivar < atmos_mean_.cols(); ivar++) {
        Array<double, 2> data_eofs(num_layers, num_eofs);
        Range packed_r(ipacked, ipacked + num_layers - 1);
        data_eofs = eofs_(all, packed_r).transpose(secondDim, firstDim);
        result.push_back(data_eofs);
        ipacked += num_layers;
    }
    return result;
}

void PCAEigenSolverFortran::solve(const blitz::Array<double, 3>& gridded_data, int num_eofs)
{
    Range all = Range::all();

    if(gridded_data.rows() != 2) {
        throw Exception("First dimension of gridded_data must equal 2");
    }

    int nlayers = gridded_data.cols();
    int npoints = gridded_data.depth();
    int nlayers2 = nlayers*2;

    Array<double, 2> taudp(nlayers, npoints, blitz::ColumnMajorArray<2>());
    Array<double, 2> omega(nlayers, npoints, blitz::ColumnMajorArray<2>());

    taudp = gridded_data(0, all, all);
    omega = gridded_data(1, all, all);

    atmos_mean_.reference(Array<double, 2>(nlayers, 2, blitz::ColumnMajorArray<2>()));
    eofs_.reference(Array<double, 2>(num_eofs, nlayers2, blitz::ColumnMajorArray<2>()));
    prin_comps_.reference(Array<double, 2>(num_eofs, npoints, blitz::ColumnMajorArray<2>()));

    bool fail;
    int message_len = 100;
    int trace_len = 100;
    blitz::Array<char, 1> message(message_len);
    blitz::Array<char, 1> trace(trace_len);

    pca_eigensolver(
            &num_eofs, &npoints, &nlayers, &nlayers2, &num_eofs, &npoints, &nlayers, &nlayers2, 
            taudp.dataFirst(), omega.dataFirst(), 
            atmos_mean_.dataFirst(), eofs_.dataFirst(), prin_comps_.dataFirst(),
            &fail, &message_len, message.dataFirst(), &trace_len, trace.dataFirst());

    if (fail) {
        Exception err;
        std::string message_str = std::string(message(all).begin(), message(all).end());
        std::string trace_str = std::string(trace(all).begin(), trace(all).end());
        err << "Fortran pca_eigensolver failure: " 
            << message_str << ", " << trace_str;
        throw err;
    }
}
