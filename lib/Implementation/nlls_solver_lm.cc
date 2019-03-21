#include <nlls_solver_lm.h>
#include<algorithm>
#include<cassert>
#include<iomanip>


using namespace FullPhysics;
using namespace blitz;
using namespace std;
using namespace Eigen;


#define MAP_BRM_ERM( BMatrix ) Map< const Matrix<double, Dynamic, Dynamic, RowMajor> >(BMatrix.data(), BMatrix.rows(), BMatrix.cols())
#define MAP_BV_ECV( BVector ) Map< const VectorXd >(BVector.data(), BVector.size())
#define MAP_ECV_BV( EVector ) blitz::Array<double, 1>(EVector.data(), shape(EVector.size()), neverDeleteData)


//  Given val1 and val2, this function computes the cos (c)
//  and sin (s) for a rotation operation, that can rotate val2
//  completely onto val1, i. e.  s*val1+c*val2=0.
//
//     --     --     --    --     --   --
//     | c  -s |     | val1 |     |  r  |
//     |       |  *  |      |  =  |     |
//     | s   c |     | val2 |     |  0  |
//     --     --     --    --     --   --
//
void givens(double val1, double val2, double &c, double &s)
{
  double  tempDbl;

  if (!val2) {
    c = 1.0;
    s = 0.0;
  } else {
    if (fabs(val2) > fabs(val1)) {
      tempDbl = -val1 / val2;
      s = 1.0 / sqrt(1.0 + tempDbl * tempDbl);
      c = s * tempDbl;
    } else {
      tempDbl = -val2 / val1;
      c = 1.0 / sqrt(1.0 + tempDbl * tempDbl);
      s = c * tempDbl;
    }
  }
}


NLLSSolverLM::NLLSSolverLM(
  const boost::shared_ptr<NLLSProblem>& p, int max_cost_function_calls, 
  const NLLSSolverLM::Options& opt, double dx_tol_abs, double dx_tol_rel,
  double g_tol_abs, double g_tol_rel, bool vrbs )
    : NLLSSolver(p, max_cost_function_calls, vrbs),
      Dx_tol_abs(dx_tol_abs), Dx_tol_rel(dx_tol_rel),
      G_tol_abs(g_tol_abs), G_tol_rel(g_tol_rel),
      Opt(opt),
      CR_ratio(0), Lambda(0)
{
}



NLLSSolver::status_t NLLSSolverLM::test_dx_rel(
  const VectorXd& dx, const VectorXd& x, double dx_tol_rel )
{
  //return (dx.cwiseAbs().array() >= dx_tol_rel*(dx_tol_rel+x.cwiseAbs().array())).any()?CONTINUE:SUCCESS;

  for(VectorXd::Index i=0; i<dx.size(); i++)
    if( abs(dx(i)) >= dx_tol_rel*(dx_tol_rel + abs(x(i))) )
      return CONTINUE;

  return SUCCESS;
}



NLLSSolver::status_t NLLSSolverLM::test_dx_abs(
  const VectorXd& dx, double dx_tol_abs )
{
  //return (dx.lpNorm<Infinity>() < dx_tol_abs)?SUCCESS:CONTINUE;
  return (dx.lpNorm<1>() < dx_tol_abs)?SUCCESS:CONTINUE;
}



NLLSSolver::status_t NLLSSolverLM::test_dx(
  const Eigen::VectorXd& dx, const Eigen::VectorXd& x, double dx_tol_rel, double dx_tol_abs )
{
    // The values of the step-test tolerances determine
    // which tests are performed.  If both tolerances are
    // greater than zero, then both tests must be successful
    // to have convergence.
    //
    status_t stat_rel = test_dx_rel(dx, x, dx_tol_rel);
    status_t stat_abs = test_dx_abs(dx, dx_tol_abs);
    status_t status = CONTINUE;
    if( (dx_tol_abs>0) && (dx_tol_rel>0) )
      status = ((stat_rel==SUCCESS) && (stat_abs==SUCCESS))?SUCCESS:CONTINUE;
    else if(dx_tol_rel>0)
      status = stat_rel;
    else if(dx_tol_abs>0)
      status = stat_abs;
    return status;
}



NLLSSolver::status_t NLLSSolverLM::test_grad_rel(
  const VectorXd& g, const VectorXd& x, double cost, double g_tol_rel )
{
  return (g.cwiseProduct( x.cwiseAbs().cwiseMax(1.0) ).lpNorm<Infinity>() < abs(g_tol_rel*max(cost,1.0)))
         ?SUCCESS:CONTINUE;
}



NLLSSolver::status_t NLLSSolverLM::test_grad_abs(
  const VectorXd& g, double g_tol_abs )
{
  //return (g.lpNorm<Infinity>() < g_tol_abs)?SUCCESS:CONTINUE;
  return (g.lpNorm<1>() < g_tol_abs)?SUCCESS:CONTINUE;
}



void NLLSSolverLM::solve()
{
  W = MAP_BRM_ERM( P->jacobian() ).colwise().stableNorm().cwiseMax(Opt.min_W);

  VectorXd X( MAP_BV_ECV(P->parameters()) );

  stat = UNTRIED;
  //
  // Any code that may cause an error before trying
  // to take a step should be placed here.
  //
  stat = CONTINUE;

  do {

    // The following three lines are only for recording purpose.
    // Info at the initial guess (the starting point) is also 
    // recorded here.
    //
    record_cost_at_accepted_point(P->cost());
    record_accepted_point(P->parameters());
    record_gradient_at_accepted_point(P->gradient());

    // Print state right after recording.
    //
    if(verbose)
      print_state();

    // After iterate() returns, Dx (step) can have a size zero (no step
    // taken) for one of the following reasons:
    //
    //   iterate() returns ERROR, and no step is computed.
    //
    //   iterate() returns STALLED, and no step is computed.
    //
    //   iterate() returns SUCCESS, and the last computed step was
    //   small enough to assume convergence, but the step was rejected
    //   because it increased the value of the cost function.
    //
    // In any of the above cases just return; however, if a step is
    // taken, just break out of the loop so that the information at
    // the new point is recorded at the end of this method.
    //
    // Moreover, when iterate() is called, it will evaluate the cost
    // function of the problem as many times as needed (almost never
    // more than just a few times) until error encountered or a good
    // step is found even if the number of the cost function evaluations
    // exceeds the limit max_cost_f_calls.
    //
    stat = iterate();
    if(Dx.size() == 0)
      return;
    else if(stat == SUCCESS)
      break;

    // The values of the gradient-test tolerances determine which
    // tests are performed.  If both tolerances are greater than
    // zero, then both tests must be successful to have convergence.
    //
    // The gradient based convergence test could be done before
    // calling iterate() just because the initial guess could be a 
    // minimum.  However, the initial guess could also be a maximum,
    // where gradient is also zero or very close zero.  Therefore,
    // it is better to go through at least one iteration to see
    // whether or not the cost function can be decreased.
    //
    status_t stat_rel = test_grad_rel(MAP_BV_ECV(P->gradient()), MAP_BV_ECV(P->parameters()), P->cost(), G_tol_rel);
    status_t stat_abs = test_grad_abs(MAP_BV_ECV(P->gradient()), G_tol_abs);
    if( (G_tol_abs>0) && (G_tol_rel>0) )
      stat = ((stat_rel==SUCCESS) && (stat_abs==SUCCESS))?SUCCESS:CONTINUE;
    else if(G_tol_rel>0)
      stat = stat_rel;
    else if(G_tol_abs>0)
      stat = stat_abs;
    if(stat == SUCCESS)
      break;

    // Performs a relative and/or absolute step convergence test.
    //
    stat = test_dx(Dx, MAP_BV_ECV(P->parameters()), Dx_tol_rel, Dx_tol_abs);
    if(stat == SUCCESS)
      break;

  } while( (stat == CONTINUE)
           && (P->num_cost_evaluations() < max_cost_f_calls)
           && (P->message() == NLLSProblem::NONE) );

  // The following three lines are only for recording purpose.
  //
  record_cost_at_accepted_point(P->cost());
  record_accepted_point(P->parameters());
  record_gradient_at_accepted_point(P->gradient());

  // Print state right after recording.
  //
  if(verbose)
    print_state();
}



void NLLSSolverLM::print_state(ostream &ostr)
{
  std::ios oldState(nullptr);
  oldState.copyfmt(ostr);

  ostr << "Solver hm_lmder;  at point # " << num_accepted_steps()
       << ";  (|f(x)|^2)/2 = " << P->cost() << ";  status = " << status_str() << endl;

  (void) ostr.precision(15);

  ostr << "Where the point x is" << endl << fixed << setw(25) << P->parameters() << endl;
  ostr << "The gradient g(x) is" << endl << fixed << setw(25) << P->gradient() << endl << endl;

  ostr.copyfmt(oldState);
}



NLLSSolver::status_t NLLSSolverLM::iterate()
{

  //  Dx (step) will get populated if a SUCCESSful step is computed.
  //
  Dx.resize(0);

  //  Perform the rank reveling QR decomposition on the Jacobian.
  //
  ColPivHouseholderQR< Matrix<double, Dynamic, Dynamic, RowMajor> > j_QR( MAP_BRM_ERM(P->jacobian()) );

  //  Just for convenience and shorter expressions.
  //
  uint32_t n = P->expected_parameter_size();
  uint32_t r = j_QR.rank();

  //  In practice it is not very likely, but in theory it is possible to
  //  get a Jacobian that is zero.  Jacobian can be zero at a minimum, at
  //  a maximum or at a saddle point.  A zero Jacobian is either valid
  //  or an error due to the implementation of the problem.  Here, we cannot
  //  assume that it is an error in the implementation of some other code.
  //  However, a zero Jacobian gives no sense of direction, and the solver
  //  will not be able to compute a step, hence iterate() returns STALLED
  //  when Jacobian is zero at the starting point.
  //
  if(r <= 0)
    return STALLED;

  double res_norm_sqrd = MAP_BV_ECV(P->residual()).squaredNorm();
  VectorXd QTResidual = j_QR.householderQ().transpose() * MAP_BV_ECV(P->residual());
  VectorXd perm_W = j_QR.colsPermutation().transpose()*W;
  MatrixXd Rrn = j_QR.matrixR().topRows(r).template triangularView<Upper>();
  MatrixXd TRrn = Rrn*perm_W.cwiseInverse().asDiagonal();

  //  Compute the Gauss-Newton step.  If the linear system is full
  //  rank, then the unique solution to the linear least square 
  //  problem is the Gauss-Newton step.  If the system is rank
  //  deficient, then solve the system for the unique minimum-norm
  //  scaled-step.  It is very important that in the case of rank
  //  deficiency, when there are infinite solutions for the step, we
  //  solve not for the unique smallest norm solution but to solve
  //  for the unique smallest norm scaled solution.
  //
  VectorXd step;
  VectorXd scaled_step;
  double scaled_step_norm;
  if(r == n) {
    step = j_QR.solve(-MAP_BV_ECV(P->residual()));
    scaled_step = W.cwiseProduct(step);
    scaled_step_norm = scaled_step.stableNorm();
  } else {
    HouseholderQR<MatrixXd> TRrnT_QR(TRrn.transpose());
    MatrixXd TrrT( TRrnT_QR.matrixQR().topRows(r).template triangularView<Upper>() );
    VectorXd y = TrrT.transpose().triangularView<Lower>().solve(-QTResidual.head(r));
    scaled_step_norm = y.stableNorm();
    VectorXd y_augmented(n);
    y_augmented << y, VectorXd::Zero(n-r);
    scaled_step = j_QR.colsPermutation()*(TRrnT_QR.householderQ()*y_augmented);
    step = W.cwiseInverse().cwiseProduct(scaled_step);
  }
  Lambda = 0;

  //  This is the loop that tries to find an acceptable step.  An
  //  acceptable step is one that significantly reduces the value
  //  of the cost function.  A significant reduction is a reduction
  //  that is larger than data noise, round-off and truncation errors.
  //
  bool cost_reducing_step=true;
  do {

    //  If the step is not within the trust region, then compute a smaller
    //  Levenberg-Marquardt step.  Only for the first iteration of the loop,
    //  at this point the step is the Gauss-Newton step and its corresponding
    //  Levenberg-Marquardt parameter is zero.  
    //
    if( scaled_step_norm > (1.0 + Opt.tr_rad_tol)*Opt.tr_rad ) {

      //  Compute an upper-bound for updating (increasing the value of)
      //  the Lev-Mar parameter.
      //
      double upper = (TRrn.transpose()*QTResidual.head(r)).stableNorm()/Opt.tr_rad;

      //  Compute a lower-bound for updating (increasing the value of)
      //  the Lev-Mar parameter.
      //
      double lower = 0;
      if(r == n) {
        double phi_0 = scaled_step_norm-Opt.tr_rad;
        VectorXd temp =j_QR.colsPermutation().transpose()*(W.asDiagonal()*scaled_step/scaled_step_norm);
        VectorXd z = Rrn.transpose().triangularView<Lower>().solve(temp);
        double phi_prime_0 = -scaled_step_norm * z.squaredNorm();
        lower = -phi_0/phi_prime_0;
      }
      if( !cost_reducing_step )
        lower = max(Lambda, lower);

      //  At this point  (lower >= 0) && (upper > lower)  must be true.
      assert( (0 <= lower) && (lower < upper) );

      bool step_in_trust_region_found=false;
      do {

        //  Levenberg-Marquard parameter (lambda) limit check.  If lambda is
        //  not in the open interval (lower, upper), then use the algorithm
        //  to choose a value in the open interval.
        //
        if( (Lambda <= lower) || (Lambda >= upper) )
          Lambda = max(0.001*upper, sqrt(lower*upper));

        //  Compute a smaller step for lambda > 0.  Givens rotations form a
        //  major part of this step.
        //
        VectorXd QTRes_head = QTResidual.head(n);
        VectorXd QTRes_extra = VectorXd::Zero(n);
        MatrixXd Rnn = MatrixXd::Zero(n,n);
        Rnn.topRows(r) = Rrn;
        MatrixXd diag = (sqrt(Lambda)*(j_QR.colsPermutation().transpose()*W)).asDiagonal();
        for(uint32_t i=0; i<n; i++) {
          double cs, sn;

          for(uint32_t j=i; j>0; j--) {
            givens( diag(j-1,i), diag(j,i), cs, sn);
            //
            RowVectorXd tempv = diag.block(j-1,i,1,n-i);
            diag.block(j-1,i,1,n-i) = cs*tempv - sn*diag.block(j,i,1,n-i);
            diag.block(j,i,1,n-i) = sn*tempv + cs*diag.block(j,i,1,n-i);
            //
            double temps = QTRes_extra[j-1];
            QTRes_extra[j-1] = cs*temps - sn*QTRes_extra[j];
            QTRes_extra[j] = sn*temps + cs*QTRes_extra[j];
          }

          givens( Rnn(i,i), diag(0,i), cs, sn );
          //
          RowVectorXd tempv = Rnn.block(i,i,1,n-i);
          Rnn.block(i,i,1,n-i) = cs*tempv - sn*diag.block(0,i,1,n-i);
          diag.block(0,i,1,n-i) = sn*tempv + cs*diag.block(0,i,1,n-i);
          //
          double temps = QTRes_head[i];
          QTRes_head[i] = cs*temps - sn*QTRes_extra[0];
          QTRes_extra[0] = sn*temps + cs*QTRes_extra[0];
        }
        VectorXd y = Rnn.triangularView<Upper>().solve(-QTRes_head);
        step = j_QR.colsPermutation()*y;

        //  After finding the smaller Levenberg-Marquardt step, scale the 
        //  step and find its norm.
        //
        scaled_step = W.cwiseProduct(step);
        scaled_step_norm = scaled_step.stableNorm();

        //  True if the norm of the new scaled step is almost equal to the
        //  size of the trust region radius, false otherwise.
        //
        step_in_trust_region_found = 
          (scaled_step_norm >= ((1.0 - Opt.tr_rad_tol)*Opt.tr_rad)) &&  (scaled_step_norm <= ((1.0 + Opt.tr_rad_tol)*Opt.tr_rad));

        //  If the step is still outside of the trust region, then update
        //  lower, lambda and upper for another iteration to find a smaller
        //  step.
        //
        if(!step_in_trust_region_found) {
          double phi = scaled_step_norm-Opt.tr_rad;
          VectorXd temp =j_QR.colsPermutation().transpose()*(W.asDiagonal()*scaled_step/scaled_step_norm);
          VectorXd z = Rnn.transpose().triangularView<Lower>().solve(temp);
          double phi_prime = -scaled_step_norm * z.squaredNorm();

          //  Update the upper limit of the Levenberg-Marquardt parameter
          //
          if(phi < 0) 
            upper = Lambda;

          //  Update the lower limit of the Levenberg-Marquardt parameter.
          //
          lower = max(lower, (Lambda-phi/phi_prime));

          //  Update the Levenberg-Marquardt parameter.
          //
          Lambda = Lambda-((phi+Opt.tr_rad)/Opt.tr_rad)*(phi/phi_prime);
        }

      } while(!step_in_trust_region_found);

    }

    double jacob_step_norm_sqrd = 
      (Rrn.triangularView<Upper>()*(j_QR.colsPermutation().transpose()*step)).squaredNorm();

    //  Here, we save the current state of the problem before moving
    //  the problem to the next point.  In this way if the step is not
    //  acceptable, then we can easily restore the state of the problem
    //  before taking the step.
    //
    blitz::Array<double, 1> x = P->parameters();

    //  Go to the next point for trial.
    //
    blitz::Array<double, 1> x_next(x+MAP_ECV_BV(step));
    P->parameters( blitz::Array<double,1>(x_next));
    double next_res_norm_sqrd = MAP_BV_ECV(P->residual()).squaredNorm();

    //  The updating of the trust region radius depends on
    //  next_res_norm_sqrd whether or not the step is acceptable.
    //  Therefore, if the problem encounters an error, it will 
    //  not be possible to continue solving the problem.  If 
    //  error, go back to the previous point and return ERROR.
    //
    if(P->message() == NLLSProblem::ERROR) {
      P->parameters(x);
      return ERROR;
    }

    //  Compute the ratio of the actual to the predicted reduction
    //  in the value of the cost function with the current step.
    //
    if( (next_res_norm_sqrd >= res_norm_sqrd) || (scaled_step_norm <= 0) ) {
      CR_ratio = 0;
    } else {
      CR_ratio = (1.0 - next_res_norm_sqrd/res_norm_sqrd) /
        ( (jacob_step_norm_sqrd/res_norm_sqrd) 
          + (2.0*Lambda*scaled_step_norm*scaled_step_norm/res_norm_sqrd) );
    }

    //  If there is no significant reduction in the value of the cost
    //  function then the step is not good.
    //
    cost_reducing_step = (CR_ratio > Opt.cr_ratio_tol);
    if(!cost_reducing_step) {

      //  Step does not reduce the value of the cost function;
      //  therefore, go back to the previous point.
      //
      P->parameters(x);

    } else {

      //  Step is accepted; therefore, update the scaling matrix
      //  using the Jacobian at the new point.
      //
      W = W.cwiseMax( MAP_BRM_ERM(P->jacobian()).colwise().stableNorm().transpose() );

      //  However, if the problem encounters an error when evaluating
      //  the Jacobian at the new point, then go back to the previous 
      //  point and return ERROR.
      //
      if(P->message() == NLLSProblem::ERROR) {
        P->parameters(x);
        return ERROR;
      }

      //  The new step is good. It reduces the values of the cost
      //  function, and there is not Jacobian evaluation error at 
      //  the new point; therefore, save the step.
      //
      Dx = step;

    }

    //  If the step gets small enough to assume convergence, then
    //  return SUCCESS regardless whether or not the step is accepted.
    //
    if(test_dx(step, MAP_BV_ECV(x_next), Dx_tol_rel, Dx_tol_abs) == SUCCESS)
      return SUCCESS;

    //  Update the trust-region radius.
    //
    double tr_scaling_factor = 1.0;
    if(CR_ratio <= 0.25) {
      if(next_res_norm_sqrd <= res_norm_sqrd)
        tr_scaling_factor = 0.5;
      else if(next_res_norm_sqrd >= 100*res_norm_sqrd)
        tr_scaling_factor = 0.1;
      else {
        double temp = -( jacob_step_norm_sqrd/res_norm_sqrd 
                         + Lambda*scaled_step_norm*scaled_step_norm/res_norm_sqrd );
        tr_scaling_factor = temp / (2.0*temp + 1.0 - next_res_norm_sqrd/res_norm_sqrd);
        if(tr_scaling_factor < 0.1)
          tr_scaling_factor = 0.1;
        else if(tr_scaling_factor > 0.5)
          tr_scaling_factor = 0.5;
      }
      Opt.tr_rad *= tr_scaling_factor;
    } else if( (CR_ratio > 0.25 && CR_ratio < 0.75 && Lambda <= 0) || (CR_ratio >= 0.75) ) {
      tr_scaling_factor = 2.0;
      Opt.tr_rad = tr_scaling_factor*scaled_step_norm;
    }

    //  If at this point the implementor of the problem believes
    //  the problem is solved, then just return SUCCESS.
    //
    if(P->message() == NLLSProblem::SOLVED)
      return SUCCESS;

  } while(!cost_reducing_step && (P->message() == NLLSProblem::NONE));

  return CONTINUE;
}
