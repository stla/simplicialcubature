{-# LANGUAGE ForeignFunctionInterface #-}
module Export
  where
import qualified Cubature as CUB
import qualified Cubature2 as CUB2
import           Foreign
import           Foreign.C


foreign export ccall test :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test rule maxevals result = do
  rule <- peek rule
  maxevals <- peek maxevals
  ([value], [errest], nevals, fl)
    <- CUB.example2 (fromIntegral maxevals) (fromIntegral rule)
  pokeArray result [realToFrac value, realToFrac errest,
                    fromIntegral nevals, (fromIntegral.fromEnum) fl]

foreign export ccall test2 :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test2 :: Ptr CInt -> Ptr CInt -> Ptr CDouble -> IO ()
test2 rule maxevals result = do
  rule <- peek rule
  maxevals <- peek maxevals
  ([value], [errest], nevals, fl)
    <- CUB2.example2 (fromIntegral maxevals) (fromIntegral rule)
  pokeArray result [realToFrac value, realToFrac errest,
                    fromIntegral nevals, (fromIntegral.fromEnum) fl]
