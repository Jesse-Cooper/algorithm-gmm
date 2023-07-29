-- Purpose  Provides a method to find the maximum likelihood estimates (MLEs)
--          for a bivariate normal mixture model (biGMM).

-- This can be extended to any dimensional multivariate normal as mentioned in
-- `Distro.hs`.


module BiGmm (biGmm, biGmmProbXandZ, gmmProbZgivenX) where


import Distro
import Matrix


-- Max iterations and the minimum threshold for the incomplete log-likelihood to
-- improve before stopping.
maxIter = 100
epsilon = 1e-5

-- Added to each estimated covariance's diagonals to prevent it from being
-- singular.
coVarReg = 1e-6


-- | Finds the maximum likelihood estimates (MLEs) for a bivariate normal
--   mixture model (biGMM).
--
--   The estimate's order may be different for different initial values.
--
--   This is a worker wrapper for `biGmmWorker`.
biGmm :: RealFloat a
      => Matrix a
      -> Vector a
      -> Matrix a
      -> Vector (Matrix a)
      -> (Vector a, Matrix a, Vector (Matrix a))
biGmm xs weights means coVars = estimates
    where
        probXandZ = biGmmProbXandZ xs weights means coVars
        logLike   = gmmLogLike probXandZ
        estimates = biGmmWorker xs weights means coVars 0 logLike


-- | Finds the current iteration of the maximum likelihood estimates (MLEs) for
--   a bivariate normal mixture model (biGMM).
--
--   Each recursion finds the next iteration of parameters and stops when the
--   max recursions (`maxIter`) or minimum threshold for the incomplete
--   log-likelihood to improve is reached.
biGmmWorker :: RealFloat a
            => Matrix a
            -> Vector a
            -> Matrix a
            -> Vector (Matrix a)
            -> Int
            -> a
            -> (Vector a, Matrix a, Vector (Matrix a))
biGmmWorker xs weights means coVars i logLike
    | i == maxIter - 1 || delta < realToFrac epsilon
        = (newWeights, newMeans, newCoVars)
    | otherwise
        = biGmmWorker xs newWeights newMeans newCoVars (i + 1) newLogLike
    where
        probXandZ    = biGmmProbXandZ xs weights means coVars
        probZgivenX  = gmmProbZgivenX probXandZ
        obsMixture   = foldRows (+) probZgivenX

        newWeights = gmmNewWeights xs probZgivenX obsMixture
        newMeans   = gmmNewMeans   xs probZgivenX obsMixture
        newCoVars  = gmmNewCoVars  xs probZgivenX obsMixture newMeans

        newProbXandZ = biGmmProbXandZ xs newWeights newMeans newCoVars
        newLogLike   = gmmLogLike newProbXandZ
        delta        = newLogLike - logLike


-- | Calculates P(X, Z) for each X and each Z.
--
--   Find the probability (PDF) of each X belonging to each mixture, then
--   multiply each mixture's probability by their weight.
--   Each row of the matrix is a mixture of all X's probability.
biGmmProbXandZ :: RealFloat a
               => Matrix a
               -> Vector a
               -> Matrix a
               -> Vector (Matrix a)
               -> Matrix a
biGmmProbXandZ xs weights means coVars = probXandZ
    where
        probX     = zipWith (biNormPdf xs) means coVars
        probXandZ = zipColsWith (*) probX weights


-- | Calculates P(Z | X) for each X and each Z.
--
--   Finds P(X) by summing P(X, Z) each Z (columns) for each X.
--   Finds P(Z | X) by dividing each Z (row) by their P(X) for each X.
gmmProbZgivenX :: RealFloat a => Matrix a -> Matrix a
gmmProbZgivenX probXandZ = probZgivenX
    where
        probX       = foldCols (+) probXandZ
        probZgivenX = zipRowsWith (/) probXandZ probX


-- | Calculates the next iteration of each GMM mixture's weight.
--
--   The update equation is just each value of `obsMixture` divided by the
--   number of observations.
gmmNewWeights :: RealFloat a
              => Matrix a
              -> Matrix a
              -> Vector a
              -> Vector a
gmmNewWeights xs probZgivenX obsMixture = newWeights
    where
        nObs       = fromIntegral $ length xs
        newWeights = map (/nObs) $ obsMixture


-- | Calculates the next iteration of each GMM mixture's mean.
--
--   The mean's update equation has been split into 3 components:
--   1. For each mixture (row), multiply X with their corresponding probability.
--          Ax
--   2. For each mixture (row), sum the vector results from `component1`.
--          sum Ax
--   3. Divide each mixture of `component2` by their corresponding `obsMixture`.
--          mean = (sum Ax) / sum A
gmmNewMeans :: RealFloat a
            => Matrix a
            -> Matrix a
            -> Vector a
            -> Matrix a
gmmNewMeans xs probZgivenX obsMixture = newMeans
    where
        component1 = map (zipColsWith (*) xs) probZgivenX
        component2 = foldRows (zipWith (+)) component1
        newMeans   = zipColsWith (/) component2 obsMixture


-- | Calculates the next iteration of each GMM mixture's covariance.
--
--   The covariance's update equation has been split into 5 components:
--   1. For each mixture, subtract its mean from each X.
--          (x - mean)
--   2. For each mixture, perform an outer product on results in `component1`.
--          (x - mean)(x - mean)^T
--   3. For each mixture and result of `component2`, multiply them by their
--      corresponding P(X | Z).
--          A(x - mean)(x - mean)^T
--   4. For each mixture, sum each matrix in `component3`.
--          sum A(x - mean)(x - mean)^T
--   5. Divide each mixture of `component5` by their corresponding `obsMixture`.
--          coVar = (sum A(x - mean)(x - mean)^T) / sum A
--   A small value is added to each covariance's diagonals to ensure they are
--   positive matrices (det > 0).
gmmNewCoVars :: RealFloat a
             => Matrix a
             -> Matrix a
             -> Vector a
             -> Matrix a
             -> Vector (Matrix a)
gmmNewCoVars xs probZgivenX obsMixture newMeans = newCoVarsReg
    where
        component1   = map (zipRowsWith (-) xs) newMeans
        component2   = mapMatrix outerProd component1
        component3   = zipMatrixWith (mapMatrixWith (*)) component2 probZgivenX
        component4   = foldRows (zipMatrixWith (+)) component3
        newCoVars    = zipWith (mapMatrixWith (/)) component4 obsMixture
        newCoVarsReg = map (mapDiagonals (+ realToFrac coVarReg)) newCoVars


-- | Calculates the incomplete log-likelihood of a GMM.
--
--   Finds P(X) by summing P(X, Z) each Z (columns) for each X.
--   Finds the incomplete log-likelihood by summing each P(Z)'s log.
gmmLogLike :: RealFloat a => Matrix a -> a
gmmLogLike probXandZ = logLike
    where
        probX   = foldCols (+) probXandZ
        logLike = sum $ map log probX
