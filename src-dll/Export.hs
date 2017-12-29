{-# LANGUAGE ForeignFunctionInterface #-}
module Export
  where
import           Cubature
import           Foreign
import           Foreign.C


foreign export ccall test :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test rule maxevals result = do
  rule <- peek rule
  maxevals <- peek maxevals
  ([value], [errest], nevals, fl)
    <- example2 (fromIntegral maxevals) (fromIntegral rule)
  pokeArray result [realToFrac value, realToFrac errest,
                    fromIntegral nevals, (fromIntegral.fromEnum) fl]
