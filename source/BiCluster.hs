-- Purpose  Provides a method of clustering points to mixtures of a bivariate
--          normal mixture model.

-- This can be extended to any dimensional multivariate normal as mentioned in
-- `Distro.hs`.


module BiCluster (biGmmClustering) where


import BiGmm
import Matrix


-- | Determines the most probable Z each X belongs to; based on a bivariate
--   normal distribution mixture's parameters.
--
--   P(Z | X) represents the probability X belongs to Z, so the highest
--   probability is the Z an X belongs to.
--
--   Each Z (row) is labelled by its index (cluster number), then the highest
--   P(Z | X) for each X is assigned as its owner.
biGmmClustering :: RealFloat a
                => Matrix a
                -> Vector a
                -> Matrix a
                -> Vector (Matrix a)
                -> Vector (a, Int)
biGmmClustering xs weights means coVars = clusters
    where
        probXandZ   = biGmmProbXandZ xs weights means coVars
        probZgivenX = gmmProbZgivenX probXandZ
        nObs        = length probZgivenX
        labelled    = zipColsWith makePair probZgivenX [0 .. nObs - 1]
        clusters    = foldCols (maxBy fst) labelled


-- | Finds max of 2 item after applying a function to both.
maxBy :: Ord b => (a -> b) -> a -> a -> a
maxBy f x y
    | f x > f y = x
    | otherwise = y


-- | Creates a tuple pair of two items.
makePair :: a -> b -> (a, b)
makePair x y = (x, y)
