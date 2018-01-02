{-# LANGUAGE DuplicateRecordFields #-}
module Cubature
  where
import Internal
import Simplex
import qualified Data.Vector.Unboxed         as UV

data Results = Results
  { values         :: [Double]
  , errorEstimates :: [Double]
  , evaluations    :: Int
  , success        :: Bool
  } deriving Show

data Result = Result
  { value         :: Double
  , errorEstimate :: Double
  , evaluations   :: Int
  , success       :: Bool
  } deriving Show

integrateOnSimplex
    :: (UVectorD -> UVectorD) -- integrand
    -> Simplices              -- domain
    -> Int                    -- number of components
    -> Int                    -- maximum number of evaluations
    -> Double                 -- desired absolute error
    -> Double                 -- desired relative error
    -> Int                    -- integration rule: 1, 2, 3 or 4
    -> IO Results             -- integral, error, evaluations, success
integrateOnSimplex f s ncomp maxevals absError relError rule = do
  let n = length (head s) - 1
  case isValidSimplices s of
    True -> do
      v <- simplicesToArray s
      (vals, errors, nevals, fl) <-
        adsimp n ncomp maxevals f absError relError rule v False
      return $ Results (UV.toList vals) (UV.toList errors) nevals (not fl)
    False -> error "invalid simplices"

integrateOnSimplex'
    :: (UVectorD -> Double)   -- integrand
    -> Simplices              -- domain
    -> Int                    -- maximum number of evaluations
    -> Double                 -- desired absolute error
    -> Double                 -- desired relative error
    -> Int                    -- integration rule: 1, 2, 3 or 4
    -> IO Result              -- integral, error, evaluations, success
integrateOnSimplex' f s maxevals absError relError rule = do
  let n = length (head s) - 1
  case isValidSimplices s of
    True -> do
      v <- simplicesToArray s
      (val, err, nevals, fl) <-
        adsimp n 1 maxevals ((UV.singleton).f) absError relError rule v False
      return $ Result (UV.head val) (UV.head err) nevals (not fl)
    False -> error "invalid simplices"

fExample :: UVectorD -> UVectorD
fExample v = let list = UV.toList v in UV.fromList [sum list, sum (map (^2) list)]

fExample' :: UVectorD -> UVectorD
fExample' v = let list = UV.toList v in UV.fromList [sum list, sum list]

example rule = integrateOnSimplex fExample [canonicalSimplex 3] 2 10000 0 1e-5 rule

example' rule = integrateOnSimplex fExample' [canonicalSimplex 3] 2 10000 0 1e-5 rule

fExample2 :: UVectorD -> Double
-- fExample2 v = UV.singleton $ sqrt((x!!3-x!!2)/(x!!1-x!!0))*exp(-(x!!1-x!!0))
fExample2 v = exp(0.5*(log (x!!3-x!!2) - log (x!!1-x!!0)) - (x!!1-x!!0))
  where y = UV.toList v
        x = map (\i -> sum $ take i y) [1..4]

example2 maxevals rule = integrateOnSimplex' fExample2 [canonicalSimplex 4] maxevals 0 1e-5 rule

fExample2' :: UVectorD -> Double
fExample2' v = sqrt((x!!3-x!!2)/(x!!1-x!!0))*exp(-(x!!1-x!!0))
  where x = UV.toList v

example2' maxevals rule = integrateOnSimplex' fExample2'
                          [[[0,0,0,0],[1,1,1,1],[0,1,1,1],[0,0,1,1],[0,0,0,1]]]
                          maxevals 0 1e-5 rule

fExample3 :: UVectorD -> Double
fExample3 v = exp (UV.sum v)

example3 maxevals rule = integrateOnSimplex' fExample3
                         [[[0,0,0],[1,1,1],[0,1,1],[0,0,1]]]
                         maxevals 0 1e-5 rule
