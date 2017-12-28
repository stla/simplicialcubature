module Common
 where
import           Data.Array.IO               (IOUArray)

type IO3dArray = IOUArray (Int,Int,Int) Double
