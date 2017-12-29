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

foreign export ccall ptr :: Ptr (Ptr ()) -> IO ()
ptr :: Ptr (Ptr ()) -> IO ()
ptr result = do
  p <- mallocBytes (2*sizeOf(undefined::Double)) :: IO (Ptr (Double))
  pokeArray p [1, 2]
  pp <- mallocBytes (2*sizeOf(undefined::Double))
  poke pp p
  poke result $ castPtr pp
