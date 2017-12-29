{-# LANGUAGE ForeignFunctionInterface #-}
module Export
  where
import Cubature
import           Foreign
import           Foreign.C


foreign export ccall test :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test rule maxevals result = do
  rule <- peek rule
  maxevals <- peek maxevals
  ([value], _, _, _) <- example2 (fromIntegral maxevals) (fromIntegral rule)
  poke result $ realToFrac value
