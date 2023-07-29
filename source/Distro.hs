-- Purpose  Provides a method to sample from and compute the PDF of a bivariate
--          normal distribution mixture.

-- Although this module is for a bivariate normal distribution, all functions
-- can be expanded to any dimensional multivariate normal, which was outside
-- this project's scope.


module Distro (biNormPdf, sampleBiNormMix) where


import Matrix
import Rand


-- | Evaluates the bivariate normal distribution PDF for each vector
--   (observation) in the matrix.
--
--   Finds the determinant and inverse of the covariance, then evaluate the
--   bivariate normal PDF as mathematically defined.
biNormPdf :: RealFloat a => Matrix a -> Vector a -> Matrix a -> Vector a
biNormPdf xs mean coVar = probs
    where
        detCoVar  = det2x2 coVar
        invCoVar  = inv2x2 coVar
        scale     = 1 / (sqrt detCoVar * 2 * pi)
        quadForms = map (quadForm invCoVar . zipWith (-) mean) xs
        probs     = map ((*scale) . exp . negate . (/2)) quadForms


-- | Samples a given number of values from a weighted mixture of bivariate
--   normals with the mixture index each was sampled from.
--
--   Each recursion randomly selects the parameters from a mixture, then samples
--   a point from that distribution.
sampleBiNormMix :: RealFloat a
                => Int
                -> Vector a
                -> Matrix a
                -> Vector (Matrix a)
                -> Integer
                -> (Vector Int, Matrix a)
sampleBiNormMix n weights means coVars s1
    | n <= 0    = ([], [])
    | otherwise = (z:zs, x:xs)
    where
        (z, (mean, coVar), s2) = sample weights (zip means coVars) s1
        (x, s3)                = sampleBiNorm mean coVar s2
        (zs, xs)               = sampleBiNormMix (n - 1) weights means coVars s3


-- | Samples a single random value from a bivariate normal distribution along
--   with a new seed number.
--
--   Combining a sequence of independent N(0, 1) random values with the Cholesky
--   of the desired covariance can simulate multivariate normals.
--
--   There are as many N(0, 1) random values as the covariance's size.
--
--   Creates the sequence (only 2 for bivariate) of N(0, 1) random values,
--   then multiplies them as a vector with the Cholesky (A*z) and adds the
--   mean to the resulting vector. (A*z + mean)
sampleBiNorm :: RealFloat a
             => Vector a
             -> Matrix a
             -> Integer
             -> (Vector a, Integer)
sampleBiNorm mean coVar s1 = (x, s3)
    where
        cholesky = cholesky2x2 coVar
        (z1, s2) = sampleNorm 0 1 s1
        (z2, s3) = sampleNorm 0 1 s2
        prod     = foldRows (+) $ zipRowsWith (*) cholesky [z1, z2]
        x        = zipWith (+) mean prod


-- | Samples a random value from the normal distribution along and a new seed
--   number.
--
--   Uses the Box–Muller transform method to sample.
sampleNorm :: RealFloat a => a -> a -> Integer -> (a, Integer)
sampleNorm mean var s1 = (x, s3)
    where
        (u1, s2) = uniform s1
        (u2, s3) = uniform s2
        r        = var * sqrt (-2 * log u1)
        x        = r * cos (2 * pi * u2) + mean
