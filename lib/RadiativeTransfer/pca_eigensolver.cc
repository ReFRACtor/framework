// PCA using Eigenproblem methods

#include "pca_eigensolver.h"
#include "fp_exception.h"

using namespace FullPhysics;
using namespace blitz;

extern "C" {
    void prepare_eigenmatrix(int *m, int *n, int *t, double *x, double* y, double *ccm);
    void pca_asymtx(double *aad, int *m, int *ia, int *ievec, int *ievec2, double *tol, double *evecd, double *evald, int *ier, double *wkd, int *message_len, char *message, bool *bad_status);
    void pca_ranker(int *n, double* arrin, int *indx);
}

PCAEigenSolver::PCAEigenSolver(std::vector<Array<double, 2> >& gridded_data, int num_eofs)
{
    solve(gridded_data, num_eofs);
}

// Adapted from the original Fortran code to be more generic

void PCAEigenSolver::solve(std::vector<Array<double, 2> >& gridded_data, int num_eofs)
{
    Range all = Range::all();
    secondIndex i2;

    int num_vars = gridded_data.size();

    // Figure out number of data values to pack from the input data
    // Expected that the second dimension of each of the variable should
    // be the number of grid points and all variables should be the same
    int num_packed = 0;
    int num_points = 0;
    for(int ivar; ivar < num_vars; ivar++) {
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
        Range packed_r(ipacked, ipacked + num_data);
        packed_data(packed_r, all) = log(packed_data(all, all));
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
    int ier = 999;
    bool asymtx_failure = false;
    
    int message_len = 100;
    blitz::Array<char, 1> asymtx_message(message_len);

    // Solve eigenproblem using PCA_ASYMTX
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
        int iord = num_packed + 1 - ipack2;
        ksq_ordered(ipack2) = ksq_abs(order(iord)-1);

        for(int ipack1 = 0; ipack1 < num_packed; ipack1++) {
            evec_2(ipack1, ipack2)= evec(ipack1, order(iord)-1);
        }
    }

    // 3. Prepare EOF and PC output 
    // =========================
    
    // Resize EOFS and Principal Components
    eofs_.resize(num_eofs, num_packed);
    prin_comps_.resize(num_eofs, num_points);

    // *** Only perform for the first few EOFs
    for(int aa = 0; aa < num_eofs; aa++) {

        // set the usable eigenvalue
        double lambda = ksq_ordered(aa);

        // Normalize everything according to SQRT(Lambda)
        double stdv = std::sqrt(lambda);

        // EOFs (Unnormalized) --> Transpose the Eigenvectors
        eofs_(aa, all) = evec_2(all, aa);

        // Project data onto E1 basis
        // -- Set the principal components (unnormalized)
        for(int w = 0; w < num_points; w++) {
            prin_comps_(aa, w) = sum(eofs_(aa, all) * packed_data(all,w));
        }

        // Final normalization of EOFs and PCs
        eofs_(aa, all) = eofs_(aa, all) * stdv;
        prin_comps_(aa, all) = prin_comps_(aa, all) / stdv;
    }

}
