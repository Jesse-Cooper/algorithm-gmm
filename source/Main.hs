-- Compile  ghc -O2 Main.hs
-- Purpose  Demonstrates the GMM algorithm for a bivariate normal distribution.

-- The covariance matrices must be positive-definite.

-- For any EM algorithm, the estimate's order may be different for different
-- initial values. The initial values in this example have been ordered such
-- that the estimates are in the same order as the true values.

-- This can be extended to any dimensional multivariate normal as mentioned in
-- `Distro.hs`.


import BiCluster
import BiGmm
import Distro
import Matrix


-- True values
trueWeights = [0.26, 0.41, 0.33]

trueMeans = [[ -1,  0.2],
             [2.2,  5.3],
             [5.8, -2.8]]

trueCoVars = [[[   2,    1],
               [   1,  1.4]],

              [[ 0.8, -1.3],
               [-1.3,  3.4]],

              [[ 1.6,  0.1],
               [ 0.1,  1.7]]]

-- Initial values
initWeights = [0.6, 0.26, 0.14]

initMeans = [[6.2,    0],
             [3.7, -2.4],
             [ 10,  6.1]]

initCoVars = [[[   3,    0],
               [   0,  0.1]],

              [[ 0.5, -2.2],
               [-2.2,   10]],

              [[ 2.7,  0.6],
               [ 0.6,  5.2]]]

-- Random seed to sample with and the number of points to sample
seed = 7907
nObs = 500

-- Number of decimal places to round results
sigFigs = 3


-- | Counts the number of correct Z (clusters) for each X.
similarClusters :: Ord a => [a] -> Vector (b, a) -> Int
similarClusters [] [] = 0
similarClusters (x:xs) ((_, y):ys)
    | x == y    = 1 + similarClusters xs ys
    | otherwise = similarClusters xs ys


-- | Shows a matrix of float numbers by first rounding the numbers.
showMatrix :: (Integral a, RealFloat b, Show b) => a -> Matrix b -> String
showMatrix _ [] = ""
showMatrix sg (row:rows) = rowStr ++ "\n" ++ showMatrix sg rows
    where
        rowStr = show $ map (roundFigs sg) row


-- | Rounds a floating point number to a given number of significant figures.
roundFigs :: (Integral a, RealFloat b) => a -> b -> b
roundFigs sg num = (fromIntegral . round $ num * multiplier) / multiplier
    where
        multiplier = 10^sg


main :: IO ()
main = do

    -- Show true values
    putStrLn "True weights:"
    print $ map (roundFigs sigFigs) trueWeights
    putStrLn "\nTrue means:"
    putStrLn $ showMatrix sigFigs trueMeans
    putStrLn "True covariances:"
    mapM_ (putStrLn . showMatrix sigFigs) trueCoVars

    -- Sample mixtures
    let (zs, xs) = sampleBiNormMix nObs trueWeights trueMeans trueCoVars seed

    -- Estimate parameters
    let (weights, means, coVars) = biGmm xs initWeights initMeans initCoVars

    -- Cluster and check the accuracy of the model
    let clusters = biGmmClustering xs weights means coVars
    let correctClusters = similarClusters zs clusters
    let accuracy = fromIntegral correctClusters / fromIntegral nObs * 100

    putStrLn $ "From " ++ (show nObs) ++ " observations"

    -- Show estimated values
    putStrLn "\nEstimated weights:"
    print $ map (roundFigs sigFigs) weights
    putStrLn "\nEstimated means:"
    putStrLn $ showMatrix sigFigs means
    putStrLn "Estimated covariances:"
    mapM_ (putStrLn . showMatrix sigFigs) coVars

    putStrLn "Clustering inference accuracy:"
    putStrLn $ show accuracy ++ "%"
