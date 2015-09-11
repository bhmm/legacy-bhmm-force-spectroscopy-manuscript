/**
 * Bayesian HMM helper functions.
 *
 * @author John D. Chodera.
 */

import java.util.Random;

public class bhmm_helper {

  /**
   * Compute log of emission probability for a normal emission model.
   *
   * @param o observation
   * @param mu mean
   * @param sigma standard deviation
   * @return the log emission probability of observation o ~ N(mu, sigma^2)
   */
  public static double logEmissionProbability(double o, double mu, double sigma) {
    double chi = ((o - mu) / sigma);
    double log_P = - (1./2.)*Math.log(2.*Math.PI) - Math.log(sigma) - (1./2.)*(chi*chi);
    return log_P;
  }

  /**
   * Determine the maximum value of a one-dimensional array.
   * 
   * @param A the array
   * @return the minimum value
   */
  static double max(double [] A) {
    double max = A[0];
    for(int i = 1; i < A.length; i++)
      if(A[i] > max)
        max = A[i];
      
    return max;
  }

  /**
   * Compute the log of a sum of exponentials in a numerically stable manner, given the arguments of the exponentials.
   * 
   * @param exp_arg_i the arguments of the exponentials to be summed.
   * @return the log of the exponential sum.
   */
  public static double logSum(double [] exp_arg_i) {
    // Determine maximum exp arg.
    double max_exp_arg = max(exp_arg_i);
      
    // Accumulate sum of exponentials after subtracting max_exp_arg to ensure the maximum argument is zero, and hence the maximum exp(x) value is unity.
    double sum = 0.0;
    for(int i = 0; i < exp_arg_i.length; i++)
      sum += Math.exp(exp_arg_i[i] - max_exp_arg);
    
    // Compute log of sum, adding back contribution from max_exp_arg.
    double logSum = Math.log(sum) + max_exp_arg;
    
    // Return log of the sum.
    return logSum;
  }  

  /**
   * Draw an integer 0...(N-1) from the given stationary distribution.
   *
   * @param p_i   p_i[i] is the probability of outcome i
   * @return an integer 0..(N-1) drawn with probability p_i[i]
   */
  public static int draw(double [] p_i, Random random) {
    // Draw an outcome.
    double r = random.nextDouble();
    int outcome;
    for(outcome = 0; outcome < p_i.length; outcome++) {
      if(r <= p_i[outcome]) break;
      r -= p_i[outcome];
    }

    // Return outcome.
    return outcome;      
  }


  /**
   * Sample a state trajectory given the observations and model.
   *
   * @param o_t   trajectory of observations
   * @param mu    mu[i] is mean of emission model of state i
   * @param sigma sigma[i] is standard deviation of emission model of state i
   * @param logPi logPi[i] is log stationary probability of state i
   * @param logTij logTij[i][j] is the log conditional transition probability from state i to state j
   * @return a sampled state assignment
   */
  public static int [] sampleStateTrajectory(double [] o_t, double [] mu, double [] sigma, double [] logPi, double [][] logTij) {

    // Initialize a new random number generator.
    Random random = new Random();
    
    // Determine timeseries length.
    int T = o_t.length;
    
    // Determine number of states.
    int nstates = mu.length;

    // Forward algorithm.

    double [][] log_alpha_it = new double[nstates][T];
    double [] exp_arg_i = new double[nstates];

    for(int i = 0; i < nstates; i++) {
      log_alpha_it[i][0] = logPi[i] + logEmissionProbability(o_t[0], mu[i], sigma[i]);
    }

    for(int t = 1; t < T; t++)
      for(int j = 0; j < nstates; j++) {        
        for(int i = 0; i < nstates; i++) 
          exp_arg_i[i] = log_alpha_it[i][t-1] + logTij[i][j] + logEmissionProbability(o_t[t], mu[j], sigma[j]);
        log_alpha_it[j][t] = logSum(exp_arg_i);
      }

    /*
    System.out.printf("log_alpha_it = \n");
    for(int t = 0; t < 10; t++) {
      for(int i = 0; i < nstates; i++)
        System.out.printf("%10.4f", log_alpha_it[i][t]);
      System.out.printf("\n");
    }
    */

    // Sample trajectory, working backwards.
    
    int [] s_t = new int[T]; // state trajectory - s_t[t] is state at time t (indexed from 1)
    
    // Draw final state.
    double [] log_p_i = new double[nstates];
    double log_normalization;
    double [] p_i = new double[nstates];

    for(int i = 0; i < nstates; i++)
      log_p_i[i] = log_alpha_it[i][T-1];
    log_normalization = logSum(log_p_i);
    for(int i = 0; i < nstates; i++)
      p_i[i] = Math.exp(log_p_i[i] - log_normalization);
    s_t[T-1] = draw(p_i, random); // indexed from 1
    
    // Work backwards
    for(int t = T-2; t >= 0; t--) {
      for(int i = 0; i < nstates; i++)
        log_p_i[i] = log_alpha_it[i][t] + logTij[i][s_t[t+1]];
      log_normalization = logSum(log_p_i);
    for(int i = 0; i < nstates; i++)
      p_i[i] = Math.exp(log_p_i[i] - log_normalization);
    
    s_t[t] = draw(p_i, random); // indexed from 1
    }

    // Increment all states by 1.
    for(int t = 0; t < T; t++)
      s_t[t] = s_t[t] + 1;

    // Return state trajectory
    return s_t;    
  }

  /**
   * Forward algorithm.
   *
   * @param o_t   trajectory of observations
   * @param mu    mu[i] is mean of emission model of state i
   * @param sigma sigma[i] is standard deviation of emission model of state i
   * @param logPi logPi[i] is log stationary probability of state i
   * @param logTij logTij[i][j] is the log conditional transition probability from state i to state j
   * @param log_rho_i log_rho_i[i] is the log of the initial probability for starting in trace i
   * @return log_alpha_ti[t][i]
   */
  public static double [][] forwardAlgorithm(double [] o_t, double [] mu, double [] sigma, double [] logPi, double [][] logTij, double [] log_rho_i) {

    // Determine timeseries length.
    int T = o_t.length;
    
    // Determine number of states.
    int nstates = mu.length;

    // Scratch storage.
    double [] exp_arg_i = new double[nstates];

    // Forward algorithm.
    double [][] log_alpha_ti = new double[T][nstates];

    for(int i = 0; i < nstates; i++) 
      log_alpha_ti[0][i] = log_rho_i[i] + logEmissionProbability(o_t[0], mu[i], sigma[i]);

    for(int t = 1; t < T; t++)
      for(int j = 0; j < nstates; j++) {        
        for(int i = 0; i < nstates; i++) 
          exp_arg_i[i] = log_alpha_ti[t-1][i] + logTij[i][j];
        log_alpha_ti[t][j] = logSum(exp_arg_i) + logEmissionProbability(o_t[t], mu[j], sigma[j]);
      }
    
    // Return alpha.
    return log_alpha_ti;
  }

  /**
   * Backward algorithm.
   *
   * @param o_t   trajectory of observations
   * @param mu    mu[i] is mean of emission model of state i
   * @param sigma sigma[i] is standard deviation of emission model of state i
   * @param logPi logPi[i] is log stationary probability of state i
   * @param logTij logTij[i][j] is the log conditional transition probability from state i to state j
   * @param log_rho_i log_rho_i[i] is the log of the initial probability for starting in trace i
   * @return log_beta_ti[t][i]
   */
  public static double [][] backwardAlgorithm(double [] o_t, double [] mu, double [] sigma, double [] logPi, double [][] logTij, double [] log_rho_i) {

    // Determine timeseries length.
    int T = o_t.length;
    
    // Determine number of states.
    int nstates = mu.length;

    // Scratch storage.
    double [] exp_arg_i = new double[nstates];

    // Backward algorithm.
    double [][] log_beta_ti = new double[T][nstates];

    for(int i = 0; i < nstates; i++) {
      log_beta_ti[T-1][i] = 0.0;
    }

    for(int t = T-2; t >= 0; t--)
      for(int i = 0; i < nstates; i++) {        
        for(int j = 0; j < nstates; j++) 
          exp_arg_i[j] = logTij[i][j] + logEmissionProbability(o_t[t+1], mu[j], sigma[j]) + log_beta_ti[t+1][j];
        log_beta_ti[t][i] = logSum(exp_arg_i);
      }

    // Return beta.
    return log_beta_ti;    
  }

  /**
   * Baum-Welch computation of transition expectations.
   *
   * @param o_t   trajectory of observations
   * @param mu    mu[i] is mean of emission model of state i
   * @param sigma sigma[i] is standard deviation of emission model of state i
   * @param logPi logPi[i] is log stationary probability of state i
   * @param logTij logTij[i][j] is the log conditional transition probability from state i to state j
   * @param log_alpha_ti[t][i] is log of alpha computed by forward algorithm
   * @param log_beta_ti[t][i] is log of beta computed by backward algorithm
   * @return log_xi_tij[t][i][j] is the probability a transition from state i to state j occurred at time t
   */
  public static double [][][] baumWelch_xi(double [] o_t, double [] mu, double [] sigma, double [] logPi, double [][] logTij, double [][] log_alpha_ti, double [][] log_beta_ti) {

    // Determine timeseries length.
    int T = o_t.length;
    
    // Determine number of states.
    int nstates = mu.length;

    // Compute log-probability of observation sequence o_t given model (used subsequently as a normalizing constant)
    double log_O = logSum(log_alpha_ti[T-1]);

    // Compute log transition probabilities using Baum-Welch.
    double [][][] log_xi_tij = new double[T-1][nstates][nstates];
    for(int t = 0; t < T-1; t++) 
      for(int i = 0; i < nstates; i++)
        for(int j = 0; j < nstates; j++)
          log_xi_tij[t][i][j] = log_alpha_ti[t][i] + logTij[i][j] + logEmissionProbability(o_t[t+1], mu[j], sigma[j]) + log_beta_ti[t+1][j] - log_O;
    
    // Return log transition expectations.
    return log_xi_tij;    
  }

  /**
   * Baum-Welch computation of log symbol output probabilities.
   *
   * @param o_t   trajectory of observations
   * @param mu    mu[i] is mean of emission model of state i
   * @param sigma sigma[i] is standard deviation of emission model of state i
   * @param logPi logPi[i] is log stationary probability of state i
   * @param logTij logTij[i][j] is the log conditional transition probability from state i to state j
   * @param log_alpha_ti[t][i] is log of alpha computed by forward algorithm
   * @param log_beta_ti[t][i] is log of beta computed by backward algorithm
   * @return log_gamma_ti[t][i] is the log probability symbol o_t[t] was emitted from state i at time t
   */
  public static double [][] baumWelch_gamma(double [] o_t, double [] mu, double [] sigma, double [] logPi, double [][] logTij, double [][] log_alpha_ti, double [][] log_beta_ti) {

    // Determine timeseries length.
    int T = o_t.length;
    
    // Determine number of states.
    int nstates = mu.length;

    // Scratch storage.
    double [] exp_arg_i = new double[nstates];

    // Compute log transition probabilities using Baum-Welch.
    double [][] log_gamma_ti = new double[T][nstates];
    for(int t = 0; t < T; t++)
      for(int i = 0; i < nstates; i++) {
        for(int j = 0; j < nstates; j++) 
          exp_arg_i[j] = log_alpha_ti[t][j] + log_beta_ti[t][j];       
        double denom = logSum(exp_arg_i);
        log_gamma_ti[t][i] = log_alpha_ti[t][i] + log_beta_ti[t][i] - denom;
      }
    
    // Return log transition expectations.
    return log_gamma_ti;    
  }

  /**
   * Calculation of emission probabilities for each state.
   *
   * @param o_t   trajectory of observations
   * @param mu    mu[i] is mean of emission model of state i
   * @param sigma sigma[i] is standard deviation of emission model of state i
   * @param Pi    Pi[i] is the normalized weight or equilibrium probability of state i
   * @returns p_ti p_it[t][i] is the probability that o_t[t] came from state i
   */
  public static double [][] computeStateProbabilities(double [] o_t, double [] mu, double [] sigma, double [] Pi) {

    // Determine timeseries length.
    int T = o_t.length;
    
    // Determine number of states.
    int nstates = mu.length;

    // Determine log weights of states.
    double [] log_Pi = new double[nstates];
    for(int i = 0; i < nstates; i++)
      log_Pi[i] = Math.log(Pi[i]);
    
    // Allocate storage for p_it
    double [][] log_p_ti = new double[T][nstates];

    // Compute log emission probabilities.
    for(int t = 0; t < T; t++) {    
      for(int i = 0; i < nstates; i++)
        log_p_ti[t][i] = log_Pi[i] + logEmissionProbability(o_t[t], mu[i], sigma[i]);
      double log_denom = logSum(log_p_ti[t]);
      for(int i = 0; i < nstates; i++)
        log_p_ti[t][i] -= log_denom;
    }

    // Compute p_it.
    double [][] p_ti = new double[T][nstates];
    for(int t = 0; t < T; t++) {    
      for(int i = 0; i < nstates; i++)
        p_ti[t][i] = Math.exp(log_p_ti[t][i]);
      // Normalize.
      double denom = 0.0;
      for(int i = 0; i < nstates; i++)      
        denom += p_ti[t][i];
      for(int i = 0; i < nstates; i++)      
        p_ti[t][i] /= denom;
    }

    return p_ti;
  }

  /**
   * Compute log-likelihood of observed state sequence given data, assuming that the first state is given.
   *
   * @param s_t   trajectory of hidden states
   * @param o_t   trajectory of observations
   * @param mu    mu[i] is mean of emission model of state i
   * @param sigma sigma[i] is standard deviation of emission model of state i
   * @param logPi logPi[i] is log stationary probability of state i
   * @param logTij logTij[i][j] is the log conditional transition probability from state i to state j
   * @param equilibrium   if True, assumes initial sample is from equilibrium
   * @return the log-likelihood of the observed state sequence
   */
  public static double computeTrajectoryLogLikelihood(int [] s_t, double [] o_t, double [] mu, double [] sigma, double [] logPi, double [][] logTij, boolean equilibrium) {
    // Determine timeseries length.
    int T = o_t.length;    
    
    // Determine number of states.
    int nstates = mu.length;

    double log_likelihood = 0.0;

    // Initialize log-probability with emission probability from first state.
    {
      int i = s_t[0] - 1;
      if (equilibrium) log_likelihood += logPi[i];      
      log_likelihood += logEmissionProbability(o_t[0], mu[i], sigma[i]);
    }

    // Compute the probability of this sequence of transitions.
    for(int t = 1; t < T; t++) {
      int i = s_t[t-1] - 1; // previous state
      int j = s_t[t] - 1; // current state
      log_likelihood += logTij[i][j] + logEmissionProbability(o_t[t], mu[j], sigma[j]);
    }

    // Return log likelihood.
    return log_likelihood;
  }

  /**
   * Compute log-likelihood of observed data given model.
   *
   * @param o_t   trajectory of observations
   * @param mu    mu[i] is mean of emission model of state i
   * @param sigma sigma[i] is standard deviation of emission model of state i
   * @param logPi logPi[i] is log stationary probability of state i
   * @param logTij logTij[i][j] is the log conditional transition probability from state i to state j
   * @param equilibrium   if True, assumes initial sample is from equilibrium
   * @return the log-likelihood of the observed data given the model
   */
  public static double computeModelLogLikelihood(double [] o_t, double [] mu, double [] sigma, double [] logPi, double [][] logTij, boolean equilibrium) {
    // Determine timeseries length.
    int T = o_t.length;    
    
    // Determine number of states.
    int nstates = mu.length;

    double log_likelihood = 0.0;

    // Compute log_rho_i.
    double [] log_rho_i = new double[nstates];
    for (int i = 0; i < nstates; i++) {
      if (equilibrium)
        log_rho_i[i] = logPi[i];
      else
        log_rho_i[i] = - Math.log((float)nstates);
    }

    // Compute forward algorithm.
    double [][] log_alpha_ti = forwardAlgorithm(o_t, mu, sigma, logPi, logTij, log_rho_i);

    // Compute total likelihood.
    log_likelihood = logSum(log_alpha_ti[T-1]);

    // Return log likelihood.
    return log_likelihood;
  }


/**
 * Compute log-likelihood of gaussian mixture model
 *
 * @param o_t   trajectory of observations
 * @param mu    mu[i] is mean of emission model of state i
 * @param sigma sigma[i] is standard deviation of emission model of state i
 * @param logPi logPi[i] is log stationary probability of state i
 * @return the log-likelihood of the observations
 */
public static double computeGmmLogLikelihood(double [] o_t, double [] mu, double [] sigma, double [] logPi) {

	
	// Determine timeseries length.
int T = o_t.length;    

// Determine number of states.
int nstates = mu.length;
	double [] exp_arg_i = new double[nstates];

double log_likelihood = 0.0;

// Compute the probability of sequence of observations
for(int t = 1; t < T; t++) {
  for(int i = 0; i < nstates; i++) 
    exp_arg_i[i] = logPi[i] + logEmissionProbability(o_t[t], mu[i], sigma[i]);
  log_likelihood += logSum(exp_arg_i);
}

// Return log likelihood.
return log_likelihood;
}

}