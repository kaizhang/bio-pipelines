{-# LANGUAGE FlexibleContexts       #-}
{-# LANGUAGE FlexibleInstances      #-}
{-# LANGUAGE FunctionalDependencies #-}
{-# LANGUAGE MultiParamTypeClasses  #-}
{-# LANGUAGE OverloadedStrings      #-}
{-# LANGUAGE TemplateHaskell        #-}

module Bio.Pipeline.CallPeaks
    ( CallPeakOpts(..)
    , CallPeakOptSetter
    , tmpDir
    , qValue
    , gSize
    , pair
    , defaultCallPeakOpts
    , callPeaks
    ) where

import           Bio.Data.Experiment.Types
import           Control.Lens
import           Control.Monad.State.Lazy
import qualified Data.Text                 as T
import           Shelly                    hiding (FilePath)
import           System.IO.Temp            (withTempDirectory)

import           Bio.Pipeline.Utils

type CallPeakOptSetter = State CallPeakOpts ()

data CallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir :: FilePath
    , callPeakOptsQValue :: Double
    , callPeakOptsGSize  :: String
    , callPeakOptsPair   :: Bool
    --, callPeakOptsBroad :: Bool
    --, callPeakOptsBroadCutoff :: Double
    }

makeFields ''CallPeakOpts

defaultCallPeakOpts :: CallPeakOpts
defaultCallPeakOpts = CallPeakOpts
    { callPeakOptsTmpDir = "./"
    , callPeakOptsQValue = 0.01
    , callPeakOptsGSize = "mm"
    , callPeakOptsPair  = False
    --, callPeakOptsBroad = False
    --, callPeakOptsBroadCutoff = 0.05
    }

-- | Call peaks using MACS2.
callPeaks :: FilePath           -- ^ Ouptut file
          -> FileSet            -- ^ Sample
          -> Maybe FileSet      -- ^ Input/control sample
          -> CallPeakOptSetter  -- ^ Options
          -> IO FileSet
callPeaks output target input setter = do
    macs2 output (target^._Single.location) (fmap (^._Single.location) input)
        fileFormat opt
    return $ Single $ format .~ NarrowPeakFile $ location .~ output $
        tags .~ ["macs2"] $ emptyFile
  where
    opt = execState setter defaultCallPeakOpts
    fileFormat
        | opt^.pair = case target of
            Single fl -> case fl^.format of
                BedFile -> "BEDPE"
                BedGZip -> "BEDPE"
                _ -> error "Only BED input is supported in pairedend mode."
            _ -> error "Paired files detected."
        | otherwise = "AUTO"
{-# INLINE callPeaks #-}

macs2 :: FilePath        -- ^ Output
      -> FilePath        -- ^ Target
      -> Maybe FilePath  -- ^ Input
      -> String          -- ^ File format
      -> CallPeakOpts
      -> IO ()
macs2 output target input fileformat opt = withTempDirectory (opt^.tmpDir)
    "tmp_macs2_dir." $ \tmp -> shelly $ do
        run_ "macs2" $
            [ "callpeak", "-f", T.pack fileformat, "-g", T.pack $ opt^.gSize
            , "--outdir", T.pack tmp, "--tempdir", T.pack tmp, "--keep-dup"
            , "all", "-q", T.pack $ show $ opt^.qValue, "-t", T.pack target
            ] ++ control
        mv (fromText $ T.pack $ tmp ++ "/NA_peaks.narrowPeak") $ fromText $
            T.pack output
  where
    control = case input of
        Nothing -> []
        Just x -> ["-c", T.pack x]
{-# INLINE macs2 #-}
