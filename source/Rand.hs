-- Purpose  Provides methods to simulate randomness.


module Rand (uniform, sample) where


-- same parameter values as C
modulus    = 2^31
multiplier = 1_103_515_245
increment  = 12_345


-- | Generates a random uniform int [0, `modulus`].
--
--   Uses a linear congruential generator.
random :: Integer -> Integer
random seed = (multiplier * seed + increment) `mod` modulus


-- | Generates a random decimal (0, 1) along with a new seed number.
--
--   `random` generates an int [0, `modulus`], adding 1 and then dividing by
--   `modulus` + 1 normalises the number to be (0, 1).
uniform :: RealFloat a => Integer -> (a, Integer)
uniform s1 = (num, s2)
    where
        s2  = random s1
        num = fromIntegral (s2 + 1) / fromIntegral (modulus + 1)


-- | Samples a weighted element from a list along with its index and a new seed
--   number.
--
--   Linearly traverses the list with a cumulative sum of the weights until a
--   uniform number is less than one of the accumulates.
sample :: RealFloat a => [a] -> [b] -> Integer -> (Int, b, Integer)
sample weights xs s1
    | sum weights < 1 = error "weights must sum to 1"
    | otherwise       = (i, x, s2)
    where
        (num, s2)  = uniform s1
        weightSums = scanl1 (+) weights
        (_, x, i)  = findBy ((<) num . fst3) (zip3 weightSums xs [0..])


-- | Finds the first element in a list that satisfies the condition.
--
--   Linearly traverses the list until the first element that satisfies the
--   condition is found.
findBy :: (a -> Bool) -> [a] -> a
findBy _ [] = error "not found"
findBy f (x:xs)
    | f x       = x
    | otherwise = findBy f xs


-- | Gets the first element of a 3 length tuple.
fst3 :: (a, b, c) -> a
fst3 (x, _, _) = x
