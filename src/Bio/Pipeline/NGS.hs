{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Bio.Pipeline.NGS where

import           Bio.Data.Experiment.Types
import           Control.Lens
import           Control.Monad             (forM)
import           Control.Monad.State.Lazy
import qualified Data.Text                 as T
import           Turtle                    hiding (FilePath, format)
import qualified Turtle                    as T

data BWAOpts = BWAOpts
    { _bwaCores    :: Int          -- ^ number of cpu cores
    , _bwaSeedLen  :: Int     -- ^ seed length, equivalent to -l
    , _bwaMaxMis   :: Int    -- ^ max mismatches in seed, equivalent to -k
    , _bwaReadTrim :: Int       -- dynamic read trimming, equivalent to -q
    } deriving (Show)

makeLenses ''BWAOpts

defaultBWAOpts :: BWAOpts
defaultBWAOpts = BWAOpts
    { _bwaCores = 1
    , _bwaSeedLen = 32
    , _bwaMaxMis = 2
    , _bwaReadTrim = 5
    }

type BWAOptSetter = State BWAOpts ()

-- | Tag alignment with BWA aligner.
bwaAlign :: FilePath  -- ^ directory to save the results
         -> FilePath  -- ^ genome index
         -> BWAOptSetter
         -> [Experiment]
         -> IO [Experiment]
bwaAlign dir' index' setter = mapM $ \e -> do
    mktree dir
    let fls = filter (\x -> x^.format == FastqFile || x^.format == FastqGZip) $ e^.files
    newFiles <- forM fls $ \fl -> do
        let output = fromText $ T.format (fp%"/"%s%"_rep"%d%".bam") dir
                (e^.eid) (fl^.replication)
            input = fromText $ T.pack $ fl^.location

        stats <- with (mktempfile "./" "bwa_align_tmp_file_XXXXXX.sai") $ \tmp -> do
            shells ( T.format (
                "bwa aln -q "%d%" -l "%d%" -k "%d%" -t "%d%" "%fp%" "%fp%" > "%fp )
                (opt^.bwaReadTrim) (opt^.bwaSeedLen) (opt^.bwaMaxMis) (opt^.bwaCores)
                index input tmp ) empty
            shells ( T.format (
                "bwa samse "%fp%" "%fp%" "%fp%
                " | samtools view -Su - | samtools sort - "%fp )
                index tmp input output ) empty
            mv (fromText $ T.format (fp%".bam") output) output
            (i, stats) <- shellStrict (T.format ("samtools flagstat "%fp) output) empty
            case i of
                ExitSuccess -> return stats
                _ -> error "samtools error!"

        return $ info .~ [("stat", stats)] $
                 format .~ BamFile $
                 location .~ T.unpack (T.format fp output) $
                 keywords .~ ["raw bam file"] $ fl

    return $ files .~ newFiles $ e
  where
    opt = execState setter defaultBWAOpts
    index = fromText $ T.pack index'
    dir = fromText $ T.pack dir'


-- | Remove low quality and redundant tags.
filterBam :: FilePath  -- ^ directory to save the results
          -> [Experiment] -> IO [Experiment]
filterBam dir' = mapM $ \e -> do
    mktree dir
    let fls = filter (\x -> x^.format == BamFile) $ e^.files
    newFiles <- forM fls $ \fl -> do
        let output = fromText $ T.format (fp%"/"%s%"_rep"%d%".bam") dir
                (e^.eid) (fl^.replication)
            input = fromText $ T.pack $ fl^.location

        shells ( T.format ("samtools view -F 0x70c -q 30 -b "%s%" > "%fp)
            input output ) empty

        return $ info .~ [] $
                 format .~ BamFile $
                 location .~ T.unpack (T.format fp output) $
                 keywords .~ ["filtered bam file"] $ fl
    return $ files .~ newFiles $ e
  where
    dir = fromText $ T.pack dir'
