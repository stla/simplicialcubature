module SphericalSimplexCubature
  (integrateOnSphericalSimplex, SphericalSimplex, orthants, Result (..))
  where
import           Simplex
import           SimplexCubature
import           SphericalSimplexCubature.Internal

integrateOnSphericalSimplex
    :: ([Double] -> Double)   -- integrand
    -> SphericalSimplex       -- domain
    -> Int                    -- maximum number of evaluations
    -> Double                 -- desired absolute error
    -> Double                 -- desired relative error
    -> Int                    -- integration rule: 1, 2, 3 or 4
    -> IO Result              -- integral, error, evaluations, success
integrateOnSphericalSimplex f ssimplex = integrateOnSimplex' f' [simplex]
  where
    f' = transformedIntegrand ssimplex f
    simplex = canonicalSimplex (length ssimplex - 1)
