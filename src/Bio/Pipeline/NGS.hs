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
        let output = fromText $ T.format (fp%"/"%s%"_rep"%d%"_filt.bam") dir
                (e^.eid) (fl^.replication)
            input = fromText $ T.pack $ fl^.location

        shells ( T.format ("samtools view -F 0x70c -q 30 -b "%fp%" > "%fp)
            input output ) empty

        return $ info .~ [] $
                 format .~ BamFile $
                 location .~ T.unpack (T.format fp output) $
                 keywords .~ ["filtered bam file"] $ fl
    return $ files .~ newFiles $ e
  where
    dir = fromText $ T.pack dir'

-- | Remove duplicates
removeDuplicates :: FilePath -> FilePath -> [Experiment] -> IO [Experiment]
removeDuplicates picardPath dir' = mapM $ \e -> do
    mktree dir
    let fls = filter (\x -> x^.format == BamFile) $ e^.files
    newFiles <- forM fls $ \fl -> do
        let output = fromText $ T.format (fp%"/"%s%"_rep"%d%"_filt_mono.bam") dir
                (e^.eid) (fl^.replication)
            input = fromText $ T.pack $ fl^.location
            qcFile = fromText $ T.format (fp%"/"%s%"_rep"%d%"_picard.qc") dir
                (e^.eid) (fl^.replication)

        with (mktempfile "./" "picard_tmp_file_XXXXXX.bam") $ \tmp -> do
            -- Mark duplicates
            shells ( T.format ("java -Xmx4G -jar "%s%" MarkDuplicates INPUT="%fp%
                " OUTPUT="%fp%" METRICS_FILE="%fp%
                " VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false")
                (T.pack picardPath) input tmp qcFile) empty

            -- Remove duplicates. Index final position sorted BAM
            shells (T.format ("samtools view -F 0x70c -b "%fp%" > "%fp) tmp output) empty
            shells (T.format ("samtools index "%fp) output) empty
            (_, stats) <- shellStrict (T.format ("samtools flagstat "%fp) output) empty

            -- Compute library complexity.
            -- Sort by position and strand.
            -- Obtain unique count statistics
            -- PBC File output
            -- TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab]
            -- TwoReadPairs [tab] NRF=Distinct/Total [tab]
            -- PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
            (_, qc) <- shellStrict (T.format ("bedtools bamtobed -i "%fp%
                " | awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf \"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\", mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'") tmp
                ) empty

            let finalBam = format .~ BamFile $
                           info .~ [("stat", stats), ("PBC_QC", qc)] $
                           keywords .~ ["processed bam file"] $
                           location .~ T.unpack (T.format fp output) $ fl
                finalBai = format .~ BaiFile $
                           info .~ [] $
                           keywords .~ ["processed bam index file"] $
                           location .~ T.unpack (T.format (fp%".bai") output) $ fl
                dupQC = format .~ Other $
                        info .~ [] $
                        keywords .~ ["picard qc file"] $
                        location .~ T.unpack (T.format fp qcFile) $ fl

            return [finalBam, finalBai, dupQC]

    return $ files .~ concat newFiles $ e
  where
    dir = fromText $ T.pack dir'
