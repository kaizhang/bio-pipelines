{-# LANGUAGE OverloadedLists   #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE TemplateHaskell   #-}

module Bio.Pipeline.NGS
    ( BWAOpts
    , BWAOptSetter
    , bwaCores
    , bwaSeedLen
    , bwaMaxMis
    , bwaReadTrim
    , defaultBWAOpts
    , bwaMkIndex
    , bwaAlign
    , filterBam
    , removeDuplicates
    , bamToBed

    , STAROpts
    , STAROptSetter
    , starMkIndex
    , starAlign
    ) where

import           Bio.Data.Experiment.Types
import           Control.Lens
import           Control.Monad             (forM)
import           Control.Monad.State.Lazy
import qualified Data.Text                 as T
import           Turtle                    hiding (FilePath, format, stderr)
import qualified Turtle                    as T
import System.IO (hPutStrLn, stderr)


--------------------------------------------------------------------------------
-- DNA-seq
--------------------------------------------------------------------------------

data BWAOpts = BWAOpts
    { _bwaCores    :: Int          -- ^ number of cpu cores
    , _bwaSeedLen  :: Int     -- ^ seed length, equivalent to -l
    , _bwaMaxMis   :: Int    -- ^ max mismatches in seed, equivalent to -k
    , _bwaReadTrim :: Int       -- ^ dynamic read trimming, equivalent to -q
    , _bwaTmpDir   :: FilePath  -- ^ temp dir
    } deriving (Show)

makeLenses ''BWAOpts

defaultBWAOpts :: BWAOpts
defaultBWAOpts = BWAOpts
    { _bwaCores = 1
    , _bwaSeedLen = 32
    , _bwaMaxMis = 2
    , _bwaReadTrim = 5
    , _bwaTmpDir = "./"
    }

type BWAOptSetter = State BWAOpts ()

-- | Generate BWA genome index
bwaMkIndex :: FilePath
           -> FilePath   -- ^ Index prefix, e.g., /path/genome.fa
           -> IO FilePath
bwaMkIndex input prefix = do
    fileExist <- testfile (fromText $ T.pack prefix)
    if fileExist
        then hPutStrLn stderr "BWA index exists. Skipped."
        else do
            mktree $ fromText dir
            shells (T.format ("ln -s "%s%" "%s) (T.pack input) (T.pack prefix)) empty
            hPutStrLn stderr "Generating BWA index"
            shells (T.format ("bwa index -p "%s%" -a bwtsw "%s)
                (T.pack prefix) (T.pack input)) empty
    return prefix
  where
    dir = fst $ T.breakOnEnd "/" $ T.pack prefix

-- | Tag alignment with BWA aligner.
bwaAlign :: IsDNASeq a
         => FilePath  -- ^ directory to save the results
         -> FilePath  -- ^ genome index
         -> BWAOptSetter
         -> [Experiment a]
         -> IO [Experiment a]
bwaAlign dir' index' setter = mapM $ \e -> do
    mktree dir
    let fls = filter (\x -> x^.format == FastqFile || x^.format == FastqGZip) $ e^.files
    newFiles <- forM fls $ \fl -> do
        let output = fromText $ T.format (fp%"/"%s%"_rep"%d%".bam")
                dir (e^.eid) (fl^.replication)
            input = fromText $ T.pack $ fl^.location

        stats <- with (mktempfile (fromText $ T.pack $ opt^.bwaTmpDir) "bwa_align_tmp_file_XXXXXX.sai") $ \tmp -> do
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
filterBam :: IsDNASeq a
          => FilePath  -- ^ directory to save the results
          -> [Experiment a] -> IO [Experiment a]
filterBam dir' = mapM $ \e -> do
    mktree dir
    let fls = filter (\x -> x^.format == BamFile) $ e^.files
    newFiles <- forM fls $ \fl -> do
        let output = fromText $ T.format (fp%"/"%s%"_rep"%d%"_filt.bam")
                dir (e^.eid) (fl^.replication)
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
removeDuplicates :: IsDNASeq a
                 => FilePath -> FilePath -> [Experiment a] -> IO [Experiment a]
removeDuplicates picardPath dir' = mapM $ \e -> do
    mktree dir
    let fls = filter (\x -> x^.format == BamFile) $ e^.files
    newFiles <- forM fls $ \fl -> do
        let output = fromText $ T.format (fp%"/"%s%"_rep"%d%"_filt_mono.bam")
                dir (e^.eid) (fl^.replication)
            input = fromText $ T.pack $ fl^.location
            qcFile = fromText $ T.format (fp%"/"%s%"_picard.qc") dir (e^.eid)

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

bamToBed :: IsDNASeq a
         => FilePath -> [Experiment a] -> IO [Experiment a]
bamToBed dir' = mapM $ \e -> do
    mktree dir
    let fls = filter (\x -> x^.format == BamFile) $ e^.files
    newFiles <- forM fls $ \fl -> do
        let output = T.format (fp%"/"%s%"bed.gz") dir $ fst $ T.breakOnEnd "." $
                snd $ T.breakOnEnd "/" $ T.pack $ fl^.location
            bedFile = format .~ BedGZip $
                      location .~ T.unpack output $ fl
        shells (T.format ("bedtools bamtobed -i "%s%" | gzip -c > "%s)
            (T.pack $ fl^.location) output) empty
        return bedFile
    return $ files .~ newFiles $ e
  where
    dir = fromText $ T.pack dir'


--------------------------------------------------------------------------------
-- RNA-seq
--------------------------------------------------------------------------------

data STAROpts = STAROpts
    { _starCmd :: T.Text
    , _starCores :: Int
    , _starTmpDir :: FilePath
    }

makeLenses ''STAROpts

type STAROptSetter = State STAROpts ()

defaultSTAROpts :: STAROpts
defaultSTAROpts = STAROpts
    { _starCmd = "STAR"
    , _starCores = 1
    , _starTmpDir = "./"
    }


-- | Create index files for STAR
starMkIndex :: FilePath   -- ^ STAR command path
            -> FilePath   -- ^ Directory used to store genome indices
            -> [FilePath] -- ^ Fastq files
            -> FilePath   -- ^ Annotation file
            -> Int        -- ^ The length of the genomic sequence
                          -- around the annotated junction to be used in
                          -- constructing the splice junctions database. Set it
                          -- to "ReadLength-1" or 100 for general purpose.
            -> IO FilePath
starMkIndex star dir fstqs anno r = do
    dirExist <- testdir (fromText $ T.pack dir)
    if dirExist
        then hPutStrLn stderr "STAR index directory exists. Skipped."
        else do
            mktree $ fromText $ T.pack dir
            hPutStrLn stderr "Generating STAR indices"
            shells cmd empty
    return dir
  where
    cmd = T.format (s%" --runThreadN 1 --runMode genomeGenerate --genomeDir "%s%
        " --genomeFastaFiles "%s%" --sjdbGTFfile "%s%" --sjdbOverhang "%d)
        (T.pack star) (T.pack dir) fstqs' (T.pack anno) r
    fstqs' = T.intercalate " " $ map T.pack fstqs

-- | Align RNA-seq raw reads with STAR
starAlign :: FilePath                    -- ^ Output directory
          -> FilePath                    -- ^ STAR genome index
          -> STAROptSetter               -- ^ Options
          -> [Experiment RNA_Seq]
          -> IO [Experiment RNA_Seq]
starAlign dir' index' setter = mapM $ \e -> do
    mktree dir
    let fls = filter (\x -> x^.format == FastqFile || x^.format == FastqGZip) $ e^.files
    newFiles <- forM fls $ \fl -> with (mktempdir
        (fromText $ T.pack $ opt^.starTmpDir) "STAR_align_tmp_dir_XXXXXX") $
        \tmp_dir -> do
            let outputGenome = fromText $ T.format (fp%"/"%s%"_rep"%d%"_genome.bam")
                    dir (e^.eid) (fl^.replication)
                outputAnno = fromText $ T.format (fp%"/"%s%"_rep"%d%"_anno.bam")
                    dir (e^.eid) (fl^.replication)
                input = fromText $ T.pack $ fl^.location

            shells ( T.format (
                s%" --genomeDir "%fp%" --readFilesIn "%fp%
                " --outFileNamePrefix "%fp%" --runThreadN "%d%
                (if fl^.format == FastqGZip then " --readFilesCommand zcat" else "")%
                " --genomeLoad NoSharedMemory --limitBAMsortRAM 0"%

                " --outFilterType BySJout"%     -- reduces the number of ”spurious” junctions
                " --outFilterMultimapNmax 20"%  -- max number of multiple alignments allowed for a read: if exceeded, the read is considered unmapped
                " --alignSJoverhangMin 8"%      -- minimum overhang for unannotated junctions
                " --alignSJDBoverhangMin 1"%    -- minimum overhang for annotated junctions
                " --outFilterMismatchNmax 999"% -- maximum number of mismatches per pair, large number switches off this filter
                " --outFilterMismatchNoverReadLmax 0.04"% -- max number of mismatches per pair relative to read length: for 2x100b, max number of mismatches is 0.06*200=8 for the paired read
                " --alignIntronMin 20"%         -- minimum intron length
                " --alignIntronMax 1000000"%    -- maximum intron length
                " --alignMatesGapMax 1000000"%  -- maximum genomic distance between mates

                " --outSAMunmapped Within  --outSAMattributes NH HI AS NM MD"%
                " --outSAMheaderCommentFile COfile.txt"%
                " --outSAMheaderHD @HD VN:1.4 SO:coordinate"%
                " --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate"%

                " --quantMode TranscriptomeSAM --sjdbScore 1"
                ) (opt^.starCmd) index input tmp_dir (opt^.starCores) ) empty

            mv (fromText $ T.format (fp%"/Aligned.sortedByCoord.out.bam") tmp_dir)
                outputGenome

            -- Sorting annotation bam
            shells ( T.format ("samtools sort -@ "%d%" -T "%fp%
                "/temp.nnnnnn.bam -o "%fp%
                " "%fp%"/Aligned.toTranscriptome.out.bam")
                (opt^.starCores) tmp_dir outputAnno tmp_dir ) empty

            -- Collect bam flagstats
            (_, stats1) <- shellStrict
                (T.format ("samtools flagstat "%fp) outputGenome) empty
            (_, stats2) <- shellStrict
                (T.format ("samtools flagstat "%fp) outputAnno) empty

            let genomeAlignFile = info .~ [("stat", stats1)] $
                    format .~ BamFile $
                    location .~ T.unpack (T.format fp outputGenome) $
                    keywords .~ ["raw bam file"] $ fl
                annoFile = info .~ [("stat", stats2)] $
                    format .~ BamFile $
                    location .~ T.unpack (T.format fp outputAnno) $
                    keywords .~ ["raw anno file"] $ fl
            return [genomeAlignFile, annoFile]

    return $ files .~ concat newFiles $ e
  where
    opt = execState setter defaultSTAROpts
    dir = fromText $ T.pack dir'
    index = fromText $ T.pack index'
